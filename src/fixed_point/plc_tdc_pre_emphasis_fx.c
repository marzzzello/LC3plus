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


void processPreEmphasis_fx(Word32 *d2_fx, Word16 *d2_fx_exp, Word16 fs_idx, Word16 n_bands, Word16 frame_dms, Word8 *scratchBuffer)
{
    Word16        s;
    Word32        nrg;
    Word16        smax;
    Counter       band;
    const Word16 *pre_emph;
    const Word16 *pre_emph_e;
    Word16 *      d2_band_fx_exp;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processPreEmphasis_fx", sizeof(struct {
                   Word16        s;
                   Word32        nrg;
                   Word16        smax;
                   Counter       band;
                   const Word16 *pre_emph;
                   const Word16 *pre_emph_e;
                   Word16 *      d2_band_fx_exp;
               }));
#endif

    d2_band_fx_exp = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_BANDS_NUMBER_PLC = 160 bytes */

    pre_emph   = lpc_lin_pre_emphasis[fs_idx];
    pre_emph_e = lpc_lin_pre_emphasis_e[fs_idx];
    SWITCH (frame_dms)
    {
        case 25: 
            pre_emph   = lpc_lin_pre_emphasis_2_5ms[fs_idx];
            pre_emph_e = lpc_lin_pre_emphasis_e_2_5ms[fs_idx];
            BREAK;
        case 50:
            pre_emph   = lpc_lin_pre_emphasis_5ms[fs_idx];
            pre_emph_e = lpc_lin_pre_emphasis_e_5ms[fs_idx];
            BREAK;
    }
	
	ASSERT(n_bands==20 || n_bands==40 || n_bands==60 || n_bands ==80);

    /* start processing */
    smax = -31; move16();

    FOR (band = 0; band < n_bands; band++)
    {
        nrg = Mpy_32_16(d2_fx[band], pre_emph[band]);

        if (nrg == 0)
        {
            s = 31; move16();
        }

        if (nrg != 0)
        {
            s = norm_l(nrg);
        }

        d2_fx[band]          = L_shl_pos(nrg, s);        move32();
        d2_band_fx_exp[band] = sub(pre_emph_e[band], s); move16();

        smax = s_max(smax, d2_band_fx_exp[band]);
    }

/* Rescale band energies */
    FOR (band = 0; band < n_bands; band++)
    {
        d2_fx[band] = L_shr_pos(d2_fx[band], s_min(sub(smax, d2_band_fx_exp[band]), 31)); move32();
    }
    /* Save common exponent for all bands */
    *d2_fx_exp = add(*d2_fx_exp, smax); move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


