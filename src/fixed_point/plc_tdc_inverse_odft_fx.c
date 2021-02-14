/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"

#include "functions.h"


void processInverseODFT_fx(Word32 *r_fx, Word16 *r_fx_exp, Word32 *d2_fx, Word16 d2_fx_exp, Word16 n_bands,
                           Word16 lpc_order, Word8 *scratchBuffer)
{
    Counter       i;
    Word16        s;
    Word16        n_bands2;
    Word32 *      x;
    const Word32 *inv_odft_twiddle_re;
    const Word32 *inv_odft_twiddle_im;
    Word8 *       buffer_BASOP_rfftN;


#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processInverseODFT_fx", sizeof(struct {
                   Counter       i;
                   Word16        s;
                   Word16        n_bands2;
                   Word32 *      x;
                   const Word32 *inv_odft_twiddle_re;
                   const Word32 *inv_odft_twiddle_im;
                   Word8 *       buffer_BASOP_rfftN;
                   Word32 *      params[2];
               }));
#endif

    x                  = scratchAlign(scratchBuffer, 0);                     /* Size = 320 bytes */
    buffer_BASOP_rfftN = scratchAlign(x, sizeof(*x) * (MAX_BANDS_NUMBER_PLC + MAX_BANDS_NUMBER_PLC/2)); /* Size = 480 bytes */

    ASSERT(lpc_order <= M);
    ASSERT(n_bands == 80 || n_bands == 60 || n_bands == 40 || n_bands == 20);

    n_bands2 = shr_pos_pos(n_bands, 1);

    test();
    IF (sub(n_bands, 20) == 0 || sub(n_bands, 60) == 0)
    {
      /* sort input samples */
        FOR (i = 0; i < n_bands2; i++)
        {
            x[2*i]         = d2_fx[2 * i];                 move32();
            x[2*i+1]       = 0;                            move32();
            x[n_bands + 2*i] = d2_fx[n_bands - 1 - 2 * i]; move32();
            x[n_bands + 2*i + 1] = 0;                      move32();
        }
        BASOP_cfft(&x[0], &x[1], n_bands, 2, &d2_fx_exp, (Word32*)buffer_BASOP_rfftN);
    }
    ELSE
    {
        /* sort input samples */
        FOR (i = 0; i < n_bands2; i++)
        {
            x[i]            = d2_fx[2 * i];               move32();
            x[n_bands2 + i] = d2_fx[n_bands - 1 - 2 * i]; move32();
        }

        BASOP_rfftN(x, n_bands, &d2_fx_exp, buffer_BASOP_rfftN);
    }

    inv_odft_twiddle_re = inv_odft_twiddle_80_re;
    inv_odft_twiddle_im = inv_odft_twiddle_80_im;
    IF (sub(n_bands, 20) == 0)
    {
        inv_odft_twiddle_re = inv_odft_twiddle_20_re;
        inv_odft_twiddle_im = inv_odft_twiddle_20_im;
    }
    ELSE IF (sub(n_bands, 40) == 0)
    {
        inv_odft_twiddle_re = inv_odft_twiddle_40_re;
        inv_odft_twiddle_im = inv_odft_twiddle_40_im;
    }
    ELSE IF (sub(n_bands, 60) == 0)
    {
        inv_odft_twiddle_re = inv_odft_twiddle_60_re;
        inv_odft_twiddle_im = inv_odft_twiddle_60_im;
    }

    s = norm_l(x[0]);

    /* imag[0] is always zero */
    r_fx[0] = L_shl_pos(x[0], s); move32();

    /* r_fx[0] = r_fx[0] * 1.0001 */
    r_fx[0] = Mpy_32_32(r_fx[0], 0x4001A36E); move32();
    IF (norm_l(r_fx[0]) > 0)
    {
        r_fx[0] = L_shl_pos(r_fx[0], 1);
    }
    ELSE
    {
        s = sub(s, 1);
    }

    /* post-twiddle */
    FOR (i = 1; i <= lpc_order; i++)
    {
        r_fx[i] = L_add(Mpy_32_32(L_shl(x[2 * i], s), inv_odft_twiddle_re[i - 1]),
                        Mpy_32_32(L_shl(x[2 * i + 1], s), inv_odft_twiddle_im[i - 1])); move32();
    }

    *r_fx_exp = sub(d2_fx_exp, s); move16();

    /* r_fx[0] must not be zero */
    IF (r_fx[0] == 0)
    {
        r_fx[0] = (Word32)0x7FFFFFFF; move32();
        FOR (i = 1; i <= lpc_order; i++)
        {
            r_fx[i] = 0; move32();
        }
        *r_fx_exp = 0; move16();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


