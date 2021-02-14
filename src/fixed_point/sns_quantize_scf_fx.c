/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"



static Word16 stage1_base(                    /* o  :  idx                                 */
                          const Word16 *t,    /* i  :  target SCFs                         */
                          const Word16 *cdbk, /* i  :  SCFs cdbk                           */
                          const Word16  R     /* i  :  number of rows in codebook       */
)
{
    Counter row;
    Word16  k_ptr, idx;
    Word32  L_min_mse, L_mse;
    Counter col;
    Word16  err;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("stage1_base", sizeof(struct {
                   Counter row, col;
                   Word16  k_ptr, idx, err;
                   Word32  L_min_mse, L_mse;
               }));
#endif
    BASOP_sub_sub_start("stage1_base");

/* find first vector error energy for  */
/* loop probably works with saturation , but it should not occur anyway */
    L_min_mse = L_add(0, 0);          /*  init acc with absolute  min mse sofar */
    FOR (col = 0; col < M / 2; col++) /* fixed to 8 elements */
    {
        err = sub(cdbk[col], t[col]); /* cdbk max abs value is 2048 = 2.^11 , max nb col is 2^3  max target is approx
                                         similar (2.^14/M)*2  = +/- 2048 , errmax is 4096   */
        L_min_mse = L_mac0(L_min_mse, err, err); /*  max L_min_mse is 8*4096*4096 =2.^(3+12+12) =  2.^27  */
    }

    idx = 0; move16();

    k_ptr = M / 2; move16(); /* ptr init to second row */
    FOR (row = 1; row < R; row++)
    {
        /* loop probably works with saturation , but it should not occur anyway */

        L_mse = L_add(L_min_mse, 0);      /* init acc with min mse sofar , */
        FOR (col = 0; col < M / 2; col++) /* fixed to 8 elements */
        {
            err   = sub(cdbk[k_ptr++], t[col]);
            L_mse = L_msu0(L_mse, err,
                           err); /* NB subtraction  from best MSE error  sofar in acc , saturation may not occur */
        }

        L_min_mse = L_sub(L_min_mse, L_max(L_mse, 0L)); /* ALWAYS update best MSE  error sofar    */

        if (L_mse > 0L) /*  if acc value  still is positive a new lower error energy vector was found in this row   */
        {
            idx = row; move16(); /* update  1-8 bits idx  */
        }

        /* this inner loop(always updating L_min_mse),          (L_msu, if )    consumes AV 19, WC  ~20 STL  cycles  ,
                                      compared to a conventional(L_mac, IF( ) )          AV 21  WC  ~23 STL  cycles per
           loop  */
    }
    ASSERT(idx >= 0 && idx < R);

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif

    return idx;
}

static void first_stage_split_search(const Word16 *cbk_LF, const Word16 *cbk_HF, const Word16 *target,
                                     const Word16 nbCbkEntries, Word16 *idxLF, Word16 *idxHF)
{
    /* find  base index for  each   SCF split  */
    *idxLF = stage1_base(target, cbk_LF, nbCbkEntries);
    *idxHF = stage1_base((&target[M / 2]), cbk_HF, nbCbkEntries);
}

static void processDeQuantize_stage1ScfDecStage1_fx(const Word16 *cbk_LF, const Word16 *cbk_HF, Word16 st1_idx0,
                                                    Word16 st1_idx1, Word16 *st1_vector)
{
    Counter col;
    Word16 offset0, offset1;
#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processDeQuantize_stage1ScfDecStage1_fx", sizeof(struct {
                   Word16  offset0, offset1;
                   Counter col;
               }));
#endif

    offset0 = shl_pos(st1_idx0, 3); /* mult by M/2 */
    offset1 = shl_pos(st1_idx1, 3);
    FOR (col = 0; col < M / 2; col++)
    {
        st1_vector[col]     = cbk_LF[offset0++]; move16();
        st1_vector[col + 8] = cbk_HF[offset1++]; move16();
    }
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static void processQuantize_stage1ScfEncStage1_fx(const Word16 *target_st1, Word16 *st1_vector, Word16 *st1_idx0Ptr,
                                                  Word16 *st1_idx1Ptr)

{
    BASOP_sub_sub_start("processQuantize_stage1ScfEncStage1_fx");

    first_stage_split_search(st1SCF0_7_base5_32x8_Q14, st1SCF8_15_base5_32x8_Q14, target_st1, SCF_STAGE1_NBCDKENTRIES,
                             st1_idx0Ptr, st1_idx1Ptr);

    processDeQuantize_stage1ScfDecStage1_fx(st1SCF0_7_base5_32x8_Q14, st1SCF8_15_base5_32x8_Q14, *st1_idx0Ptr,
                                            *st1_idx1Ptr, st1_vector);

    BASOP_sub_sub_end();

    return;
}

/* gain-shape MSE search in warped SCF-residual domain,  synthesis in SCF resiudal domain allows for easy weighting */

static void pvq_enc_find_best_submode_pre_post_fx(
    const Word16 *target_st2, /* this target is in the linearized  warped domain , same as input to PVQ search  */
    const Word16 *enc_pulses_far, Word16 *enc_pulses_near, const Word16 *enc_pulsesA, const Word16 *enc_pulsesB,
    Word16 *sub_mode_ptr, Word16 *i_gain_ptr, Word16 *enc_adj_glob_warped_vec, Word8 *scratchBuffer) /* Size = 18 * M */
{

    Counter       L_section, idx;
    const Word16 *search_en1shape[N_SCF_SHAPES_ST2];
    const Word16 *search_gainTab[N_SCF_SHAPES_ST2];
    Word16        search_n_gains[N_SCF_SHAPES_ST2];
    Word32        L_mse, L_mse_min, L_idx;
    Word16 *      pulses_far, *pulses_near, *pulsesAB, *pulsesA;
    Word16 *      target_w, *shape_far, *shape_near, *shapeAB, *shapeA;
    Word16  tmp, err;
    Counter i;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("pvq_enc_find_best_submode_pre_post_fx", sizeof(struct {
                   Counter i, L_section, idx;
                   Word16 *search_en1shape[N_SCF_SHAPES_ST2];
                   Word16 *search_gainTab[N_SCF_SHAPES_ST2];
                   Word16  search_n_gains[N_SCF_SHAPES_ST2];
                   Word32  L_mse, L_mse_min, L_idx;
                   Word16 *pulses_far, *pulses_near, *pulsesAB, *pulsesA;
                   Word16 *target_w, *shape_far, *shape_near, *shapeAB, *shapeA;
                   Word16  tmp, err;
               }));
#endif

    pulses_near = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * M */

    pulsesAB = (Word16 *)scratchAlign(pulses_near, sizeof(*pulses_near) * M); /* Size = 2 * M */

    pulsesA = (Word16 *)scratchAlign(pulsesAB, sizeof(*pulsesAB) * M); /* Size = 2 * M */

    target_w = (Word16 *)scratchAlign(pulsesA, sizeof(*pulsesA) * M); /* Size = 2 * M */

    shape_near = (Word16 *)scratchAlign(target_w, sizeof(*target_w) * M); /* Size = 2 * M */

    shapeAB = (Word16 *)scratchAlign(shape_near, sizeof(*shape_near) * M); /* Size = 2 * M */

    shapeA = (Word16 *)scratchAlign(shapeAB, sizeof(*shapeAB) * M); /* Size = 2 * M */

    pulses_far = (Word16 *)scratchAlign(shapeA, sizeof(*shapeA) * M); /* Size = 2 * M */

    shape_far = (Word16 *)scratchAlign(pulses_far, sizeof(*pulses_far) * M); /* Size = 2 * M */

    BASOP_sub_sub_start("pvq_enc_find_best_submode_pre_post_fx");

    /* construct pulse vectors and en1 normalized shape vectors  */ /* use shape Q in Q14 */
    basop_memmove(pulses_far, enc_pulses_far, M * sizeof(*pulses_far));
    basop_memmove(pulses_near, enc_pulses_near, M * sizeof(*pulses_near));
    basop_memmove(target_w, target_st2, M * sizeof(*target_w));

    pvq_dec_en1_normQ14_fx(shape_near, pulses_near, sns_Kval[2][0], M); /* near outlier mode  */
    pvq_dec_en1_normQ14_fx(shape_far, pulses_far, sns_Kval[3][0], M);   /* far outlier mode  */

    /* regular mode(with a split),   prepare vectors  of full length M */
    basop_memmove(pulsesAB, enc_pulsesA, N_SETA * sizeof(*pulsesAB));
    basop_memmove(pulsesA, enc_pulsesA, N_SETA * sizeof(*pulsesA));

    FOR (i = N_SETA; i < M; i++)
    {
        pulsesAB[i] = enc_pulsesB[sub(i, N_SETA)]; move16();
    }

    IF (M > N_SETA)
    {
        basop_memset(&pulsesA[N_SETA], 0, (M - N_SETA) * sizeof(*pulsesA));
    }

    pvq_dec_en1_normQ14_fx(shapeAB, pulsesAB, sns_Kval[0][0], M);
    /* regular AB , b_pulses = 1 ;*/ /* OPT: combine  with shapeA */

    pvq_dec_en1_normQ14_fx(shapeA, pulsesA, sns_Kval[1][0], M);
    /* regular A ,  b_pulses = 0 */ /* OPT:  M-> N_SETA */

    /* setup search structure */

    /* now aligned with order of  j  {regular=0, regular_lf=1, outlier_near=2, outlier far=3}  */

    search_en1shape[0] = shapeAB;
    search_gainTab[0]  = sns_gaintabPtr[0];
    search_n_gains[0]  = sns_gainSz[0]; /* assumes whole bits */

    search_en1shape[1] = shapeA;
    search_gainTab[1]  = sns_gaintabPtr[1];
    search_n_gains[1]  = sns_gainSz[1]; /* assumes whole bits */

    search_en1shape[2] = shape_near;
    search_gainTab[2]  = sns_gaintabPtr[2];
    search_n_gains[2]  = sns_gainSz[2]; /*assume whole bits */

    search_en1shape[3] = shape_far;
    search_gainTab[3]  = sns_gaintabPtr[3];
    search_n_gains[3]  = sns_gainSz[3]; /*assume whole bits */

    /* start actual search loop */

    /* basic raw MSE loop,   */
    L_mse_min = INT_MAX;         move32();
    L_idx     = L_deposit_l(-1); /* section in low 2  bits* gain idx above */

    FOR (L_section = 0; L_section < N_SCF_SHAPES_ST2; L_section++)
    {
        /* raw MSE  over gain and shape */
        FOR (idx = 0; idx < search_n_gains[L_section]; idx++)
        {
            /* MSE ( proc_target_local[i]-adjGain[i]*en1Shape[i] ) */

            L_mse = L_deposit_l(0);
            FOR (i = 0; i < M; i++)
            {
                tmp   = mult_r(search_gainTab[L_section][idx], search_en1shape[L_section][i]); /* Q15+14+1-16= Q14 */
                err   = sub(target_w[i], tmp);                                                 /*  both in  Q14      */
                L_mse = L_mac0(L_mse, err, err);                                               /* Q14+14 = Q28 */
            }

            IF (L_sub(L_mse, L_mse_min) < 0) /* OPT: always update L_mse_min) */
            {
                L_mse_min = L_mse;                          move32();
                L_idx     = L_mac0(L_section, idx, 1 << 2); /* save both section and gain  idx */
            }
        } /* gains */
    }     /*submodes*/

    L_section = L_and(0x3L, L_idx); /* section was stored in two lowest bits */
    ASSERT(L_section >= 0 && L_section <= 3);
    *i_gain_ptr = extract_l(L_shr_pos(L_idx, 2)); /*1,2,3 bit gain */
    ASSERT(*i_gain_ptr >= 0 && *i_gain_ptr <= 7);

    /* returns a scaled and transformed vector, ___EXACTLY__ as a decoder would scale it */
    ASSERT(enc_adj_glob_warped_vec != NULL);
    {
        /* warp/rotate search result to SCF residual domain */
        idct16_fx(search_en1shape[L_section], target_w); /* fwd synthesis  warping */
        /* actual synthesis gain scaling in SCF-residual domain, for easy weighting analysis  */
        pvq_dec_scale_vec_fx(target_w, search_gainTab[L_section][*i_gain_ptr], enc_adj_glob_warped_vec);
    }

    *sub_mode_ptr = extract_l(L_section); move16(); /* 0,1,2,3 */

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    return;
}

static void processQuantize_stage2ScfEncStage2_fx(const Word16 *target_st2, Word16 *st2_vector, Word32 *L_prm_idx,
                                                  Word16 submodes, Word8 *scratchBuffer) /* Size = 26 * M + 48 */
{                                                                                        /*func */
    Word16 *proc_target, *enc_pulses_far, *enc_pulses_near, *enc_pulsesA, *enc_pulsesB;

    Word16 *pulses_fin, *pulses_proj;
    Word32  L_tmp;

    Word8 *buffer_pvq_enc_find_best_submode_pre_post_fx;

    PvqEntry_fx enc_PVQ_OA, enc_PVQ_B;
    Word16      submode, i_gain, submodeMSB, submodeLSB;
    Word32 *    L_search_corr, *L_search_en;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processQuantize_stage2ScfEncStage2_fx", sizeof(struct {
                   Word16 *proc_target, *enc_pulses_far, *enc_pulses_near, *enc_pulsesA, *enc_pulsesB;

                   Word16 *pulses_fin, *pulses_proj;
                   Word32  L_tmp;

                   Word8 *buffer_pvq_enc_find_best_submode_pre_post_fx;

                   PvqEntry_fx enc_PVQ_OA, enc_PVQ_B;
                   Word16      submode, i_gain, submodeMSB, submodeLSB;
                   Word32 *    L_search_corr, *L_search_en;
               }));
#endif

    buffer_pvq_enc_find_best_submode_pre_post_fx = scratchAlign(scratchBuffer, 0); /* Size = 18 * M */
    proc_target =
        (Word16 *)scratchAlign(buffer_pvq_enc_find_best_submode_pre_post_fx,
                               sizeof(*buffer_pvq_enc_find_best_submode_pre_post_fx) * 18 * M); /* Size = 2 * M */

    enc_pulses_near = (Word16 *)scratchAlign(proc_target, sizeof(*proc_target) * M);         /* Size = 2 * M */
    enc_pulsesA     = (Word16 *)scratchAlign(enc_pulses_near, sizeof(*enc_pulses_near) * M); /* Size = 2 * N_SETA */
    enc_pulsesB     = (Word16 *)scratchAlign(enc_pulsesA, sizeof(*enc_pulsesA) * N_SETA);    /* Size = 2 * N_SETB */
    pulses_fin = (Word16 *)scratchAlign(enc_pulsesB, sizeof(*enc_pulsesB) * N_SETB); /* Size = 2 * N_SCF_SHAPES_ST2 */
    pulses_proj =
        (Word16 *)scratchAlign(pulses_fin, sizeof(*pulses_fin) * N_SCF_SHAPES_ST2); /* Size = 2 * N_SCF_SHAPES_ST2 */
    L_search_corr =
        (Word32 *)scratchAlign(pulses_proj, sizeof(*pulses_proj) * N_SCF_SHAPES_ST2); /* Size = 4 * N_SCF_SHAPES_ST2 */
    L_search_en    = (Word32 *)scratchAlign(L_search_corr,
                                         sizeof(*L_search_corr) * N_SCF_SHAPES_ST2); /* Size = 4 * N_SCF_SHAPES_ST2 */
    enc_pulses_far = (Word16 *)scratchAlign(L_search_en, sizeof(*L_search_en) * N_SCF_SHAPES_ST2); /* Size = 2 * M */

    BASOP_sub_sub_start("processQuantize_stage2ScfEncStage2_fx");

    /* fixed setup for a given  bitrate of 38 ,  no  moves needed */
    /* k_far  = sns_Kval[3][0]; */
    /* k_near = sns_Kval[2][0]; */
    /* kA     = sns_Kval[1][0]; */ /* regular, regular_lf */
                                   /* kB is always  1 */

    /* NB  these search indecese do not correspond exactly to specification shape_index j */

    pulses_fin[0] = sns_Kval[3][0]; /* far   6 */
    pulses_fin[1] = sns_Kval[2][0]; /* near  8 */
    pulses_fin[2] = sns_Kval[1][0]; /* section A     10 */
    pulses_fin[3] = sns_Kval[0][1]; /* section B     1 */

    pulses_proj[0] = sns_Kval[3][0];
    pulses_proj[1] = 0;
    pulses_proj[2] = 0;
    pulses_proj[3] = 0;

    /*  pre_process  */
    dct16_fx(target_st2, proc_target); /* enc analysis */

    /* get the initial four integer shape candidate vectors,  no normalization at this stage  */
    pvq_enc_search_fx(proc_target, enc_pulses_far, enc_pulses_near, enc_pulsesA, enc_pulsesB, L_search_corr,
                      L_search_en, pulses_fin, pulses_proj, M, N_SETA);

    /* scale with gains a after a  unit energy fwd transform  */
    /* apply transform to each candidate shape vector priot  to gain-shape search loop */
    submode = submodes; /* used as input solely to debug/unit test a specific shape mode  */

    /*target should be in a  linearized residual domain target */
    /* search pre, synthesis  post*/
    pvq_enc_find_best_submode_pre_post_fx(proc_target, enc_pulses_far, enc_pulses_near, enc_pulsesA, enc_pulsesB,
                                          &submode, &i_gain, st2_vector,
                                          buffer_pvq_enc_find_best_submode_pre_post_fx); /* Q14 tr out */

    /* send parameters  to multiplexor as a series/vector  of Long Words */
    /*    0 :    0..3  submode  */
    /*    1 :    0..7  gain_ind  */
    /*    2 :    0..1  LeadSign ind */
    /*    3 :    25 bit     MPVQ index    outl_near or  A  part  */
    /*    4 :    3.7 to 21 bits  MPVQ index           B  part  OR   -2  */

    L_prm_idx[0] = L_deposit_l(submode); /*  complete submode fwd'ed  to ari_codec as  0,1,2,3  */

    submodeMSB = shr_pos(submode, 1);                       /* LSB of submode , sent as main submode bit  */
    submodeLSB = s_and(submode, 0x1); /* LSB of submode  */ /*   sent via shape param  */

    /* gain, shape indicese , incl. calls to  MPVQ indexing */
    IF (submodeMSB == 0)
    { /* regular modes::   j=0(reg=AB)  or 1(reg_lf  A)  */ /* regular mode, with two or one shape indices  */

        /* assume regular_lf part ,  shape_j == 1 */
        enc_PVQ_OA =
            mpvq_index_fx(enc_pulsesA, N_SETA, sns_Kval[submode][0]); /* o : leading_sign_index, index, size, k_val */
        L_prm_idx[2] = L_deposit_l(enc_PVQ_OA.lead_sign_ind);         /*LS set A */

        ASSERT(enc_PVQ_OA.size == (UWord32)sns_MPVQ_Sz[submode][0]);
        L_prm_idx[3] = L_add(0L, (Word32)enc_PVQ_OA.index); /* MPVQ shape index set A fractional   */

        /* section B always have low indexing dynamics and is  combined into one joint single  index */
        IF (submodeLSB == 0)
        {                                                                              /* regular   AB  , shape_j == 0*/
            L_prm_idx[1] = L_deposit_l(i_gain); /* full established gain idx fwd'ed */ /*      2  possible values */
            enc_PVQ_B    = mpvq_index_fx(enc_pulsesB, N_SETB, 1);
            ASSERT(((enc_PVQ_B.size << 1)) ==
                   (sns_MPVQ_Sz[submode][1])); /*  two lowest indeces indicate all_zero B section  */

            L_tmp        = L_shl_pos((Word32)enc_PVQ_B.index, 1);            /* 2*section B  MPVQ index */
            L_prm_idx[4] = L_add(L_tmp, enc_PVQ_B.lead_sign_ind); move32(); /* add joint section B and  LS index */

            ASSERT(L_prm_idx[4] >= 0 && L_prm_idx[4] < (Word32)sns_MPVQ_Sz[submode][0]);
        }
        ELSE
        {
            L_prm_idx[1] = L_deposit_l(i_gain);
            /* MSBs of established gain idx */ /*  2 or 4   total  possible values */
            L_prm_idx[4] = L_deposit_l(-2);
        }
    }
    ELSE
    {
        /* outlier  modes   shape_j= 2(near, LSB=0) or 3(far, LSB=1)  */

        IF (submodeLSB == 0)
        {
            L_prm_idx[1] = L_deposit_l(i_gain); /* established gain idx  */ /*   4  possible values */
            enc_PVQ_OA   = mpvq_index_fx(enc_pulses_near, M,
                                       sns_Kval[submode][0]); /* o :  leading_sign_index,  index, size, k_val        */
            ASSERT(enc_PVQ_OA.size == sns_MPVQ_Sz[submode][0]);
            L_prm_idx[3] = L_add(0L, enc_PVQ_OA.index); /* MPVQ index  fractional bits */
            L_prm_idx[4] = L_deposit_l(-1);             /* no gain LSBs  */
        }
        ELSE
        {
            L_prm_idx[1] = L_deposit_l(i_gain); /* established gain idx MSBs   */ /*   all 4 or 8   possible values */
            enc_PVQ_OA   = mpvq_index_fx(enc_pulses_far, M,
                                       sns_Kval[submode][0]); /* o :  leading_sign_index,  index, size, k_val        */
            ASSERT(enc_PVQ_OA.size == sns_MPVQ_Sz[submode][0]);
            L_prm_idx[3] = L_add(0L, enc_PVQ_OA.index); /* MPVQ index  fractional bits */
            L_prm_idx[4] = L_deposit_l(-2);             /*  */
        }
        L_prm_idx[2] = L_deposit_l(enc_PVQ_OA.lead_sign_ind); /* LS shape single bit */
    }

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    return;
}

static Word16 scfdec_stage2_fx(                          /* o: ber flag */
                               const Word32 *L_prm_idx,  /* set to -1 if not used */
                               Word16 *      st2_vector, /*o: Q14 */
                               Word8 *       scratchBuffer)
{
    /*   MPVQ deindexing, gainscaling transform and transform */
    Dyn_Mem_Deluxe_In(
        Word16  submode;
        Word16  submodeLSB, submodeMSB;
        Word16  gValQ13;
        Word16  idxB;
        Word16  maxK;
        Word16  BER_dec;
        Word16 *dec_pulses;
        Word16 *dec_en1_vec;
        Word16 *dec_adj_glob_vec;
    );

    BASOP_sub_sub_start("scfdec_stage2_fx");

    dec_pulses       = (Word16 *)scratchAlign(scratchBuffer, 0);                      /* Size = 2 * M = 32 bytes */
    dec_en1_vec      = (Word16 *)scratchAlign(dec_pulses, sizeof(*dec_pulses) * M);   /* Size = 2 * M = 32 bytes */
    dec_adj_glob_vec = (Word16 *)scratchAlign(dec_en1_vec, sizeof(*dec_en1_vec) * M); /* Size = 2 * M = 32 bytes */

    /* get submode   */
    submode = extract_l(L_prm_idx[0]); /* 0..3 */

    submodeLSB = s_and(submode, 0x1);
    submodeMSB = shr_pos(submode, 1);

    /* get initial adjustment gain vector for  regular, outl_near   */
    ASSERT(L_prm_idx[1] >= 0 && L_prm_idx[1] < sns_gainSz[submode]);
    gValQ13 = sns_gaintabPtr[submode][L_prm_idx[1]];
    ASSERT(gValQ13 >= 0);

    /* gain, shape indices,  incl.calls  to MPVQ deindexing */
    IF (submodeMSB != 0)
    {
        /* outlier_near or outlier_far  mode decoding */
        maxK    = sns_Kval[submode][0]; move16();
        BER_dec = pvq_dec_deidx_fx(dec_pulses, maxK, M, extract_l(L_prm_idx[2]), (UWord32)L_prm_idx[3]);
    }
    ELSE
    { /* regular mode, with potentially two shape indices  */

        maxK    = sns_Kval[submode][0]; move16();
        BER_dec = pvq_dec_deidx_fx(dec_pulses, maxK, N_SETA, extract_l(L_prm_idx[2]), (UWord32)L_prm_idx[3]);

        IF (submodeLSB == 0)
        {
            idxB = extract_l(L_prm_idx[4]); /* 0..11 */
            ASSERT(idxB >= 0 && idxB < (Word16)sns_MPVQ_Sz[0][1]);
            BER_dec |= pvq_dec_deidx_fx(&(dec_pulses[N_SETA]), sns_Kval[submode][1], N_SETB, s_and(idxB, 0x1),
                                        (UWord32)L_deposit_l(shr_pos(idxB, 1)));
            /* maxK does not need to be increased as set B is not stacked  */
        }
        ELSE
        { /* LSB gain bit already parsed */
            ASSERT(L_prm_idx[4] < 0);
            basop_memset(&dec_pulses[N_SETA], 0, (N_SETB) * sizeof(*dec_pulses));
        }
    }

    /* normalize decoded integer vector , exactly as on encoder side !!  */
    pvq_dec_en1_normQ14_fx(dec_en1_vec, dec_pulses, maxK, M);

    idct16_fx(dec_en1_vec, dec_adj_glob_vec); /* fwd warping  in unscaled domain */

    /* scaling aligend with encoder search  */
    pvq_dec_scale_vec_fx(dec_adj_glob_vec, gValQ13, st2_vector);

    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
    return BER_dec;
}

void processSnsQuantizeScfEncoder_fx(Word16  scf[],        /* i: input scf M */
                                     Word32 *L_prm_idx,    /* o: indeces . negative == unused */
                                     Word16 *scf_q,        /* o: quantized scf M */
                                     Word8 * scratchBuffer) /* Size = 28 * M + 52 */
{
    Counter col;
    Dyn_Mem_Deluxe_In(
        Word16 *target_st2;
        Word16 *st1_idx; /* stage 1 indices */
        Word8 * buffer_processQuantize_stage2ScfEncStage2_fx;
    );

    target_st2 = (Word16 *)scratchAlign(scratchBuffer, 0);                    /* Size = 2 * M */
    st1_idx    = (Word16 *)scratchAlign(target_st2, sizeof(*target_st2) * M); /* Size = 2 * 2 */
    buffer_processQuantize_stage2ScfEncStage2_fx = (Word8 *)scratchAlign(st1_idx, sizeof(*st1_idx) * 2);
    /* Size = 26 * M + 48 */

    /* TBD needs update  */

    /* 1st stage trained VQ   */
    processQuantize_stage1ScfEncStage1_fx(scf, scf_q, &st1_idx[0], &st1_idx[1]);
    L_prm_idx[0] = L_deposit_l(st1_idx[0]);
    L_prm_idx[1] = L_deposit_l(st1_idx[1]);

/* 2nd stage PVQ-based SCF quantizer   */
    FOR (col = 0; col < M; col++)
    {
        target_st2[col] = sub(scf[col], scf_q[col]);
    }
    processQuantize_stage2ScfEncStage2_fx(target_st2, scf_q, &L_prm_idx[2], VQMODES26,   /* 0xF means all submodes */
                                          buffer_processQuantize_stage2ScfEncStage2_fx); /*  PVQ  in stage 2 */
    Dyn_Mem_Deluxe_Out();
}

Word16 processSnsQuantizeScfDecoder_fx(                                      /* o: BER flag */
                                       Word32 *L_prm_idx,                    /* i: indeces */
                                       Word16 scf_q[], Word8 *scratchBuffer) /* o:  M */
{
    Dyn_Mem_Deluxe_In(
        Word16 BER_flag;
    );

    /* Decode First Stage */
    processDeQuantize_stage1ScfDecStage1_fx(st1SCF0_7_base5_32x8_Q14, st1SCF8_15_base5_32x8_Q14,
                                            extract_l(L_prm_idx[0]), extract_l(L_prm_idx[1]), scf_q);

    /* Decode Second Stage */
    BER_flag = scfdec_stage2_fx(&(L_prm_idx[2]), scf_q, scratchBuffer);

    Dyn_Mem_Deluxe_Out();
    return BER_flag;
}

