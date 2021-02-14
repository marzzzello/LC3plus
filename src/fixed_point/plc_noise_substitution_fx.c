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


#ifndef NONBE_PLC4_ADAP_DAMP
void processPLCNoiseSubstitution_fx(Word32 spec[], Word16 spec_prev[], Word16 L_spec, Word16 nbLostFramesInRow,
                                    Word16 stabFac, Word16 frame_dms, Word16 *alpha, Word16 *cum_alpha, Word16 *seed)
#else
void processPLCNoiseSubstitution_fx(Word32 spec[], Word16 spec_prev[], Word16 L_spec)
#endif
{
#ifndef NONBE_PLC4_ADAP_DAMP
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word16  tmp16;
        Word16  spec_prev16;
    );

    SWITCH (frame_dms)
    {
    case 25: nbLostFramesInRow = shr(add(nbLostFramesInRow, 3), 2); BREAK;
    case 50: nbLostFramesInRow = shr(add(nbLostFramesInRow, 1), 1); BREAK;
    }

    /* get damping factor */
#ifndef NONBE_PLC4_BURST_TUNING
    IF (sub(nbLostFramesInRow, 4) < 0)
    {
        *alpha = add(26214 /*0.8 Q15*/, mult_r(6553 /*0.2 Q15*/, stabFac)); move16();
    }
    ELSE IF (sub(nbLostFramesInRow, 6) < 0)
    {
        *alpha = add(19660 /*0.6 Q15*/, mult_r(9830 /*0.3 Q15*/, stabFac)); move16();
    }
    ELSE IF (sub(nbLostFramesInRow, 8) < 0)
    {
        *alpha = add(16384 /*0.5 Q15*/, mult_r(13107 /*0.4 Q15*/, stabFac)); move16();
    }
    ELSE
    {
        *alpha = add(14745 /*0.45 Q15*/, mult_r(13107 /*0.4 Q15*/, stabFac)); move16();
    }
#else
    *alpha = add(26214 /*0.8 Q15*/, mult_r(6553 /*0.2 Q15*/, stabFac)); move16();
    IF (sub(nbLostFramesInRow, PLC_FADEOUT_IN_MS / 10) > 0)
    {
        *alpha = 0; move16();
    }
    ELSE IF (sub(nbLostFramesInRow, 2) > 0)
    {
        tmp16  = div_s(sub(PLC_FADEOUT_IN_MS / 10, nbLostFramesInRow), sub(PLC_FADEOUT_IN_MS / 10, sub(nbLostFramesInRow, 1)));
        *alpha = mult(*alpha, tmp16);
    }
#endif

    SWITCH (frame_dms)
    {
    case 25:
        IF (sub(*alpha, 32767) < 0)
        {
            tmp16  = 0;
            *alpha = Sqrt16(*alpha, &tmp16);  move16();
            *alpha = shl(*alpha, tmp16);
        }
		IF (sub(*alpha, 32767) < 0)
        {
            tmp16  = 0;
            *alpha = Sqrt16(*alpha, &tmp16);  move16();
            *alpha = shl(*alpha, tmp16);
        }
        BREAK;
    case 50:
        IF (sub(*alpha, 32767) < 0)
        {
            tmp16  = 0;
            *alpha = Sqrt16(*alpha, &tmp16); move16();
            *alpha = shl(*alpha, tmp16);
        }
        BREAK;
    }

    *cum_alpha = mult_r(*cum_alpha, *alpha); move16();

    tmp16 = *seed; move16();

    /* Add noise and damping */
    FOR (i = 0; i < L_spec; i++)
    {
        tmp16 = extract_l(L_mac0(16831, tmp16, 12821));

        spec_prev16 = mult(spec_prev[i], *cum_alpha);

        if (tmp16 < 0)
        {
            spec_prev16 = negate(spec_prev16);
        }

        spec[i] = L_deposit_h(spec_prev16);
    }

    /* High pass to prevent overflows */
    spec[0] = Mpy_32_16(spec[0], 6553 /* 0.2 Q15*/);  move32();
    spec[1] = Mpy_32_16(spec[1], 16384 /* 0.5 Q15*/); move32();

    *seed = tmp16; move16();
#else
    Dyn_Mem_Deluxe_In(
        Counter i;
    );

    FOR (i = 0; i < L_spec; i++)
    {
        spec[i] = L_deposit_h(spec_prev[i]);
    }

    /* High pass to prevent overflows */
    spec[0] = Mpy_32_16(spec[0], 6553 /* 0.2 Q15*/);  move32();
    spec[1] = Mpy_32_16(spec[1], 16384 /* 0.5 Q15*/); move32();
#endif /* NONBE_PLC4_ADAP_DAMP */

    Dyn_Mem_Deluxe_Out();
}



