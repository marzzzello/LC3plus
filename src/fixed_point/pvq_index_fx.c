/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"

#define SIGNBIT_FX 0x80000000u
#define SIGNBIT_SHRT_FX 0x8000


static void initOffsets_fx(Word16 dim_in, UWord32 *h_mem,
                           Word16 k_val_in) /* may be removed with tables for N=16,10,6  */
{
    UWord32 k_val_prev, k_val_curr;
    UWord32 k_val, UL_k_val_in;
#ifdef DYNMEM_COUNT
    Dyn_Mem_In("initOffsets_fx", sizeof(struct {
                   UWord32 k_val_prev, k_val_curr;
                   UWord32 k_val, UL_k_val_in;
               }));
#endif

    h_mem[0] = UL_deposit_l(0);
    h_mem[1] = UL_deposit_l(1);

    UL_k_val_in = UL_deposit_l(k_val_in);
    IF (sub(dim_in, 2) == 0)
    {
        FOR (k_val = 2; k_val <= UL_k_val_in; k_val++)
        {
            h_mem[k_val] = UL_subNsD(UL_lshl(k_val, 1), 1U); move32();
        }
        h_mem[k_val] = UL_k_val_in; move32();
    }
    ELSE
    {
        k_val_prev = UL_deposit_l(1U);
        FOR (k_val_curr = 2; k_val_curr <= UL_k_val_in; k_val_curr++)
        {
            h_mem[k_val_curr] = UL_addNsD(1U, UL_Mpy_32_32(k_val_curr, UL_lshl(k_val_prev, 1))); move32();
            k_val_prev        = UL_addNsD(k_val_curr, 0U);
        }
        h_mem[k_val_curr] = UL_Mpy_32_32(k_val_curr, k_val_prev); move32();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static void a_fwd_fx(UWord32 *a_in,   /* i/o: offsets   */
                     Word16   n_items /* i  :  items, k's  */
)
{
    UWord32  a_1, a_in0;
    Counter  i;
    UWord32 *a_in_prev_ptr;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("a_fwd_fx", sizeof(struct {
                   UWord32  a_1, a_in0;
                   Counter  i;
                   UWord32 *a_in_prev_ptr;
               }));
#endif

    a_in0 = UL_deposit_l(1);

    a_in_prev_ptr = &(a_in[-1]);
    FOR (i = 1; i <= n_items; i++)
    {
        a_1              = UL_addNsD(a_in0, UL_addNsD(a_in_prev_ptr[i], a_in[i]));
        a_in_prev_ptr[i] = a_in0; move32();
        a_in0            = UL_addNsD(a_1, 0U);
    }
    a_in_prev_ptr[i] = a_in0; move32();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static void a_bwd_fx(UWord32 *a_in,   /* i/o: offsets   */
                     Word16   n_items /* i:  n_items  */
)
{
    UWord32  a_1, a_in0;
    Counter  i;
    UWord32 *a_in_prev_ptr;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("a_bwd_fx", sizeof(struct {
                   UWord32  a_1, a_in0;
                   Counter  i;
                   UWord32 *a_in_prev_ptr;
               }));
#endif

    a_in0         = UL_deposit_l(0);
    a_in_prev_ptr = &(a_in[-1]);

    FOR (i = 1; i <= n_items; i++)
    {
        a_1              = UL_subNsD(UL_subNsD(a_in[i], a_in0), a_in_prev_ptr[i]);
        a_in_prev_ptr[i] = a_in0; move32();
        a_in0            = UL_addNsD(a_1, 0U);
    }
    a_in_prev_ptr[i] = a_in0; move32();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static void a_u_fwd_fx(UWord32 *a_u_in, Word16 k_val_in, Word16 mem_size_m1)
{
    UWord32 u_kp1_prev, u_kp1;
    UWord32 u_k_prev;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("a_u_fwd_fx", sizeof(struct {
                   UWord32 u_kp1_prev, u_kp1;
                   UWord32 u_k_prev;
               }));
#endif

    u_kp1_prev = a_u_in[mem_size_m1]; move32();
    u_k_prev   = UL_lshr(a_u_in[k_val_in], 1);

    a_fwd_fx(&a_u_in[1], k_val_in);

    u_kp1               = UL_lshr(a_u_in[k_val_in], 1);
    a_u_in[mem_size_m1] = UL_addNsD(1U, UL_addNsD(u_kp1_prev, UL_addNsD(u_k_prev, u_kp1)));

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static Word16 get_lead_sign_fx(UWord32 *ind)
{
    Word16 leading_sign;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("get_lead_sign_fx", sizeof(struct { Word16 leading_sign; }));
#endif

    leading_sign = 1; move16();
    if (UL_and(*ind, 1) != 0)
    {
        leading_sign = -1; move16();
    }
    (*ind) = UL_lshr(*ind, 1);

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif

    return leading_sign;
}

/*-------------------------------------------------------------------*
 * mind2vec_one_fx()
 *-------------------------------------------------------------------*/
static void mind2vec_one_fx(Word16  k_val_in,                           /* i:  nb unit pulses    , [ 0...K_MAX ]   */
                            Word16  leading_sign,                       /* i: leading sign  -1, 0, 1*/
                            UWord32 ind, /* i:  index                */ /* parameter could  be omitted */
                            Word16 *vec_out                             /* o:  pulse train          */
)
{
    *vec_out = (Word16)ind; /* dummy assignment to handle the common ind parameter warning  */

    if (leading_sign < 0)
    {
        k_val_in = negate(k_val_in);
    }
    *vec_out = k_val_in; move16();
}

static Word16 setval_update_sign_fx(Word16 k_delta, Word16 k_max_local, Word16 *leading_sign, UWord32 *ind_in,
                                    Word16 *vec_out)
{
    IF (k_delta != 0)
    {
        mind2vec_one_fx(k_delta, *leading_sign, *ind_in, vec_out);
        *leading_sign = get_lead_sign_fx(ind_in);
        k_max_local   = sub(k_max_local, k_delta);
    }
    return k_max_local;
}

/*-------------------------------------------------------------------*
 * mind2vec_fx()
 *-------------------------------------------------------------------*/
static void mind2vec_fx(Word16   dim_in,       /* i:  dimension        */
                        Word16   k_max_local,  /* i:  nb unit pulses   */
                        Word16   leading_sign, /* i: leading sign  */
                        UWord32  ind,          /* i:  index            */
                        Word16 * vec_out,      /* o:  pulse train      */
                        UWord32 *h_in          /* i:  offset vector   A=1+2U  */
)
{
    Counter pos;
    Word16  k_acc, k_delta;
    UWord32 UL_tmp_offset, UL_diff;
    UWord16 sgn;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("mind2vec_fx", sizeof(struct {
                   Counter pos;
                   Word16  k_acc, k_delta;
                   UWord32 UL_tmp_offset, UL_diff;
                   UWord16 sgn;
               }));
#endif

    k_acc = k_max_local; move16();
    FOR (pos = 0; pos < dim_in; pos++)
    {

        IF (ind != 0)
        {

            k_acc = k_max_local; move16();

            UL_tmp_offset = UL_addNsD(h_in[k_acc], 0U);

            UL_diff = UL_subNs(ind, UL_tmp_offset, &sgn);

            WHILE (sgn)
            {
                UL_diff = UL_subNs(ind, h_in[--k_acc], &sgn);
            }

            ind = UL_addNsD(UL_diff, 0U);

            k_delta = sub(k_max_local, k_acc);
        }
        ELSE
        {
            mind2vec_one_fx(k_max_local, leading_sign, ind, &vec_out[pos]);
            BREAK;
        }

        k_max_local = setval_update_sign_fx(k_delta, k_max_local, &leading_sign, &ind, &vec_out[pos]);

        a_bwd_fx(h_in, add(k_max_local, 1));
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

PvqEntry_fx get_size_mpvq_calc_offset_fx(                   /* o : size, dim, k_val        */
                                         Word16   dim_in,   /* i : dimension                */
                                         Word16   k_val_in, /* i : nb unit pulses           */
                                         UWord32 *h_mem     /* o : offsets                  */
)
{
    Counter     i;
    PvqEntry_fx entry;
    Word16      kp1;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("get_size_mpvq_calc_offset_fx", sizeof(struct {
                   Counter     i;
                   PvqEntry_fx entry;
                   Word16      kp1;
               }));
#endif

    entry.dim   = dim_in;   move16();
    entry.k_val = k_val_in; move16();

    entry.index         = L_deposit_l(0);
    entry.lead_sign_ind = 0; move16();

    ASSERT(dim_in <= M);
    ASSERT(tabledKMAX[dim_in] != 0);

    /* tabled values for worst case K */ /* made into table lookup for N=16, 10, 6  */
    kp1 = add(k_val_in, 1);
    FOR (i = 0; i <= kp1; i++) /* A+U copying */
    {
        h_mem[i] =
            UL_addNsD(MPVQ_offs_ptr[dim_in][i], 0U); /* a vector  copying is needed as MPVQ recursion is in place */
    }
    /* special handling of last  U offset   in k+1 column  */
    if (sub(k_val_in, tabledKMAX[dim_in]) != 0)
    {
        h_mem[kp1] = UL_lshr(h_mem[kp1], 1); /* (A+1)/2 , convert from  A(K+1) to  U(K+1)  domain */
    }
    entry.size =
        UL_addNsD(1U, UL_addNsD(h_mem[kp1], UL_lshr(h_mem[k_val_in], 1))); /* MPVQ size calc. 1 + H(K+1) + (A(K)>>1) */

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif

    return entry;
}

/*-------------------------------------------------------------------*
 * mpvq_deindex_fx()
 *-------------------------------------------------------------------*/
void mpvq_deindex_fx(                           /* o :  void                        */
                     const PvqEntry_fx *entry,  /* i :  sign_ind, index, dim, k_val */
                     UWord32 *          h_mem,  /* i :  A/U offsets                 */
                     Word16 *           vec_out /* o :  pulse train                 */
)
{
    Word16 leading_sign;
#ifdef DYNMEM_COUNT
    Dyn_Mem_In("mpvq_deindex_fx", sizeof(struct { Word16 leading_sign; }));
#endif
    BASOP_sub_sub_start("mpvq_deindex_fx");

    basop_memset(vec_out, 0, (entry->dim) * sizeof(Word16));

    leading_sign = 1; move16();
    if (entry->lead_sign_ind != 0)
    {
        leading_sign = -1; move16();
    }

    IF (entry->k_val != 0)
    {
        mind2vec_fx(entry->dim, entry->k_val, leading_sign, entry->index, vec_out, h_mem);
    }
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    BASOP_sub_sub_end();
}

/*-------------------------------------------------------------------*
 * vec2mind_two_fx()
 *-------------------------------------------------------------------*/
static void vec2mind_two_fx(const Word16 *vec_in,        /* i : PVQ  pulse train        */
                            Word16 *      k_val_out_ptr, /* o : number of unit pulses    */
                            UWord32 *     next_sign_ind, /* i/o: next sign ind           */
                            UWord32 *     ind            /* o: MPVQ index                */
)
{
    UWord32 lead_sign_ind_add;
    Word16  abs0, abs1, abs01, sptr;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("vec2mind_two_fx", sizeof(struct {
                   UWord32 lead_sign_ind_add;
                   Word16  abs0, abs1, abs01, sptr;
               }));
#endif

    abs0           = abs_s(vec_in[0]);
    abs1           = abs_s(vec_in[1]);
    abs01          = add(abs0, abs1);
    *k_val_out_ptr = abs01; move16();
    *ind           = UL_deposit_l(0);

    *next_sign_ind = UL_deposit_h(SIGNBIT_SHRT_FX);

    IF (abs01 != 0)
    {
        sptr           = 0; move16();
        *next_sign_ind = UL_deposit_l(sptr);

        test();
        IF (abs0 != 0 && abs1 != 0)
        {
            lead_sign_ind_add = UL_deposit_l(1);
            if (vec_in[1] < 0)
            {
                lead_sign_ind_add = UL_deposit_l(2);
            }
            *ind = UL_addNsD(UL_deposit_l((UWord16)lshl(sub(abs1, 1), 1)), lead_sign_ind_add);
        }
        ELSE
        {
            IF (abs0 == 0)
            {
                *ind = UL_deposit_l((UWord16)sub(lshl(abs1, 1), 1));
                sptr = 1; move16();
            }
        }

        if (vec_in[sptr] < 0)
        {
            *next_sign_ind = UL_deposit_l(1);
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static void enc_push_sign(Word16 val, UWord32 *next_sign_ind, UWord32 *index)
{
    test();
    IF ((UL_and(*next_sign_ind, SIGNBIT_FX) == 0) && (val != 0))
    {
        *index = UL_addNsD(UL_lshl(*index, 1), *next_sign_ind);
    }
    if (val < 0)
    {
        *next_sign_ind = UL_deposit_l(1);
    }
    if (val > 0)
    {
        *next_sign_ind = UL_deposit_l(0);
    }
}

static void vec2mind_fx(Word16        dim_in,        /* i :  dim                       */
                        Word16        k_val_in,      /* i :  number of unit pulses     */
                        const Word16 *vec_in,        /* i :  PVQ pulse train           */
                        UWord32 *     next_sign_ind, /* o :  pushed leading sign       */
                        UWord32 *     index,         /* o :  MPVQ index                */
                        UWord32 *     N_MPVQ_ptr,    /* o :  size(N_MPVQ(dim,K_val_in))*/
                        UWord32 *     h_mem)              /* o :  offsets                   */
{
    Counter pos;
    Word16  mem_size_m1, k_val_acc, tmp_val;
    UWord32 tmp_h;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("vec2mind_fx", sizeof(struct {
                   Counter pos;
                   Word16  mem_size_m1, k_val_acc, tmp_val;
                   UWord32 tmp_h;
               }));
#endif

    mem_size_m1    = add(k_val_in, 1);
    *next_sign_ind = UL_deposit_h(SIGNBIT_SHRT_FX);

    pos = sub(dim_in, 2);
    vec2mind_two_fx(&vec_in[pos], &k_val_acc, next_sign_ind, index);
    initOffsets_fx(3, h_mem, k_val_in);

    tmp_h = h_mem[k_val_acc]; move32();
    FOR (pos--; pos >= 0; pos--)
    {
        tmp_val = vec_in[pos]; move16();
        enc_push_sign(tmp_val, next_sign_ind, index);

        *index = UL_addNsD(*index, tmp_h);

        k_val_acc = add(k_val_acc, abs_s(tmp_val));

        IF (pos != 0)
        {
            a_u_fwd_fx(h_mem, k_val_in, mem_size_m1);
        }
        tmp_h = UL_addNsD(h_mem[k_val_acc], 0U);
    }
    *N_MPVQ_ptr = UL_addNsD(1U, UL_addNsD(UL_lshr(tmp_h, 1), h_mem[mem_size_m1])); move32();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

PvqEntry_fx mpvq_index_fx(                          /* o : leading_sign_index, index, size, k_val        */
                          const Word16 *vec_in,     /* i : signed pulse train        */
                          Word16        dim_in,     /* i : dimension                 */
                          Word16        k_val_local /* i : nb unit pulses            */
)
{
    PvqEntry_fx result;
    UWord32     h_mem[1 + KMAX_FX + 1];
    UWord32     lead_sign_ind;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("mpvq_index_fx", sizeof(struct {
                   PvqEntry_fx result;
                   UWord32     h_mem[1 + KMAX_FX + 1];
                   UWord32     lead_sign_ind;
               }));
#endif
    BASOP_sub_sub_start("mpvq_index_fx");

    ASSERT(k_val_local <= KMAX_FX);

    result.k_val = k_val_local; move16();
    result.dim   = dim_in;      move16();

    vec2mind_fx(dim_in, k_val_local, vec_in, &lead_sign_ind, &result.index, &result.size, h_mem);

    result.lead_sign_ind = u_extract_l(lead_sign_ind);

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    BASOP_sub_sub_end();

    return result;
}

