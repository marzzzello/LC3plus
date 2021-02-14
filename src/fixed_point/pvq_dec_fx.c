/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


Word16 pvq_dec_deidx_fx(                          /* out BER detected 1 , ok==0 */
                        Word16 *      y,          /* o:   decoded vector (non-scaled int)  */
                        const Word16  k_val,      /* i:   number of allocated pulses       */
                        const Word16  dim,        /* i:   Length of vector                 */
                        const Word16  LS_ind,     /* i; lS index              1 bit        */
                        const UWord32 UL_MPVQ_ind /* i; MPVQ  index                        */
)
{
    Dyn_Mem_Deluxe_In(
        Word16      BER_flag;
        UWord32     h_mem[1 + KMAX_FX + 1];
        PvqEntry_fx entry;
    );

    BER_flag = 0; move16();

    /* get_size will likely be called before this function,     as the range decoder needs the size to fetch the index
     */
    entry = get_size_mpvq_calc_offset_fx(dim, k_val, h_mem); /* TBD should be made into tables for N=16,10,6  */

    entry.lead_sign_ind = LS_ind;         move16();
    entry.index         = L_deposit_l(0); /* only  in case dim == 1 */
    IF (sub(dim, 1) != 0)
    {
        entry.index = UL_MPVQ_ind;

        /* safety check in case of bit errors */
        IF (L_sub(entry.index, entry.size) >= 0)
        {
            BER_flag    = 1;            move16();
            entry.index = 0; move16(); /* return something deterministic/valid, and LOW complex  */
        }
    }
    mpvq_deindex_fx(&entry, h_mem, y); /* actual deindexing  */

    Dyn_Mem_Deluxe_Out();
    return BER_flag;
}


void pvq_dec_scale_vec_fx(const Word16 *inQ14, Word16 adjGainQ13, Word16 *outQ14)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );
    FOR (i = 0; i < M; i++)
    {
        outQ14[i] = add(outQ14[i], mult_r(adjGainQ13, inQ14[i])); move16();
    }
    Dyn_Mem_Deluxe_Out();
}


void pvq_dec_en1_normQ14_fx(/*  Have to be used EXACTLY the same way in both  both encoder and decoder */
                            Word16 *      xq, /* o:   en1 normalized decoded vector (Q14)   */
                            const Word16 *y,  /* i:   decoded vector (non-scaled int)  */
                            const Word16  k_val_max,
                            /* i:   max possible K   in Q0 kO or kA   */ /* OPT:  not BE , use dynamic max  pulse
                                                                            amplitude */
                            const Word16 dim                             /* i:   Length of vector                 */
)
{

    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32  L_tmp;
        Word16  shift_num, shift_tot;
        Word16  isqrtQ16_local, tmp, exp, exp_shift;
        Word32  L_yy;
    );

/* energy normalization starts here */
    L_yy = L_mult0(y[0], y[0]);
    FOR (i = 1; i < dim; i++)
    {
        L_yy = L_mac0(L_yy, y[i], y[i]); /* stay in Q0 */ /* OPT: reuse some energies from PVQ linear search */
    }
    /* 16 bit */
    IF (L_sub(L_yy, SQRT_EN_MAX_FX) < 0)
    {
        ASSERT(L_yy > 4);                               /* Q16 isqrt table lookup not valid below  5  */
        isqrtQ16_local = isqrt_Q16tab[L_yy]; move16(); /* 1 cycle */
    }
    ELSE
    {
        /* about 8-9 cycles */
        exp            = 15; move16(); /* set ISqrt16() exp_in to get delta exp out near 0  when Lyy is in Q0  */
        tmp            = ISqrt16(extract_l(L_yy),
                      &exp); /* exp  out is now a delta shift with a later tmp Q15 multiplication in mind  */
        exp_shift      = add(exp, 16 - 15);   /*  up to Q16 */
        isqrtQ16_local = shl(tmp, exp_shift); /* new mantissa in a fixed  Q16  */
    }

    shift_num = norm_s(k_val_max);      /* simply account for the preknown fixed max possible pulseamp in y */
    shift_tot = sub(14 - 1, shift_num); /* upshift  to get  to Q14 */
    FOR (i = 0; i < dim; i++) /*  upshifted y[i]  used    */
    {
        L_tmp = L_mult(isqrtQ16_local, shl_pos(y[i], shift_num)); /* Q(16+0+shift_num +1    =   shift_num+1  */
        xq[i] = round_fx(L_shl(L_tmp, shift_tot)); move16();     /* Q14 ,      */
    }

    Dyn_Mem_Deluxe_Out();
}

