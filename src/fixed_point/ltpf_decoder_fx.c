/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


static void ltpf_synth_filter(Word16 *synth_ltp, Word16 *synth, Word16 length, Word16 pitch_int, Word16 pitch_fr,
                              Word16 gain, Word16 scale_fac_idx, Word16 fs_idx,
                              Word16 fade /* 0=normal, +1=fade-in, -1=fade-out */);

/*************************************************************************/


void process_ltpf_decoder_fx(Word16 *x_e, Word16 L_frame, Word16 old_x_len, Word16 fs_idx, Word16 old_y_len,
                             Word16 *old_e, Word16 *x_in, Word16 *old_x, Word16 *y_out, Word16 *old_y, Word16 ltpf,
                             Word16 ltpf_active, Word16 pitch_index, Word16 *old_pitch_int, Word16 *old_pitch_fr,
                             Word16 *old_gain, Word16 *mem_ltpf_active, Word16 scale_fac_idx, Word16 bfi,
                             Word16 concealMethod,
                             Word16 damping, Word16 *old_scale_fac_idx, Word8 *scratchBuffer)
{
    Counter i;
    Word16  gain, s, s0, s1, pitch, pitch_int, pitch_fr;
    Word16 *x, *y;
    Word16 *z;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("process_ltpf_decoder_fx", sizeof(struct {
                   Counter i;
                   Word16  gain, s, s0, s1, pitch, pitch_int, pitch_fr;
                   Word16 *x, *y;
                   Word16 *z;
               }));
#endif

    z = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = MAX_LEN / 4 + 10 */


#ifdef NONBE_PLC2_LTPF_FADEOUT_FIX
    test();
    IF ((sub(bfi, 1) == 0) && (concealMethod == 0))
#else
    test(); test();
    IF (sub(bfi, 1) == 0 && (sub(concealMethod, 2) == 0 || concealMethod == 0))
#endif
    {
        ltpf        = 0; move16();
        ltpf_active = 0; move16();

#ifndef NONBE_PLC_LTPF_FIX
        *mem_ltpf_active = 0; move16();
        bfi              = 0; move16();
#endif
    }

    /* Filter parameters */
    IF (sub(bfi, 1) != 0)
    {
        IF (ltpf == 0)
        {
            pitch_int = 0; move16();
            pitch_fr  = 0; move16();
        }
        ELSE
        {
            /* Decode pitch */
            IF (sub(pitch_index, 380) < 0)
            {
                pitch_int = shr_pos(add(pitch_index, 64), 2);
                pitch_fr  = add(sub(pitch_index, shl_pos(pitch_int, 2)), 128);
            }
            ELSE IF (sub(pitch_index, 440) < 0)
            {
                pitch_int = shr_pos(sub(pitch_index, 126), 1);
                pitch_fr  = sub(sub(shl_pos(pitch_index, 1), shl_pos(pitch_int, 2)), 252);
            }
            ELSE
            {
                pitch_int = sub(pitch_index, 283);
                pitch_fr  = 0; move16();
            }
            pitch     = add(shl_pos(pitch_int, 2), pitch_fr);
            pitch     = mult_r(shl_pos(pitch, 2), pitch_scale[fs_idx]);
            pitch_int = shr_pos(pitch, 2);
            pitch_fr  = sub(pitch, shl_pos(pitch_int, 2));
        }

        /* Decode gain */
        if (scale_fac_idx < 0)
        {
            ltpf_active = 0;
            ASSERT(!(*old_scale_fac_idx < 0 && *mem_ltpf_active == 1));
        }
        IF (ltpf_active == 0)
        {
            gain = 0; move16();
        }
        ELSE
        {
            ASSERT(scale_fac_idx >= 0);
            gain = gain_scale_fac[scale_fac_idx]; move16();
        }
    }
    ELSE
    {
        /* fix to avoid not initialized filtering for concelament 
           might be necessary in case of bit errors or rate switching */
        if (scale_fac_idx < 0) {
            if (*mem_ltpf_active && *old_scale_fac_idx>=0)
            {
                scale_fac_idx = *old_scale_fac_idx;
            }
        }
        
        ltpf_active = *mem_ltpf_active; move16();
#ifdef NONBE_PLC2_LTPF_FADEOUT_FIX
        if ((sub(concealMethod, 2) == 0))
        { /* always start fade off to save filtering WMOPS for the remaining 7.5 ms  */
            assert(bfi == 1);
            ltpf_active = 0; move16(); /*always start fade off , still maintain  *mem_ltpf_active */
        }
#endif
        pitch_int = *old_pitch_int;
        pitch_fr  = *old_pitch_fr;
        gain      = mult_r(*old_gain, damping);
    }

    test(); test();
    IF (ltpf_active == 0 && *mem_ltpf_active == 0)
    {
        /* LTPF inactive */

        basop_memmove(y_out, x_in, L_frame * sizeof(Word16));

        /* Update */
        s = sub(*old_e, *x_e);
        IF (s > 0)
        {
            basop_memmove(old_y, &old_y[L_frame], (old_y_len - L_frame) * sizeof(Word16));

            IF (sub(s, 15) > 0)
            {
                basop_memset(&old_y[old_y_len - L_frame], 0, (L_frame) * sizeof(Word16));

                basop_memset(old_x, 0, (old_x_len) * sizeof(Word16));
            }
            ELSE
            {
                FOR (i = 0; i < L_frame; i++)
                {
                    old_y[i + old_y_len - L_frame] = shr(x_in[i], s); move16();
                }
                FOR (i = 0; i < old_x_len; i++)
                {
                    old_x[i] = shr(x_in[i + L_frame - old_x_len], s); move16();
                }
            }
        }
        ELSE
        {
            IF (sub(s, -15) < 0)
            {
                basop_memset(old_y, 0, (old_y_len - L_frame) * sizeof(Word16));
            }
            ELSE
            {
                FOR (i = 0; i < old_y_len - L_frame; i++)
                {
                    old_y[i] = shl(old_y[i + L_frame], s); move16();
                }
            }

            basop_memmove(&old_y[old_y_len - L_frame], x_in, (L_frame) * sizeof(Word16));

            basop_memmove(old_x, &x_in[L_frame - old_x_len], (old_x_len) * sizeof(Word16));

            *old_e = *x_e; move16();
        }
        *old_pitch_int   = pitch_int; move16();
        *old_pitch_fr    = pitch_fr;  move16();
        *old_gain        = 0;         move16();
        *mem_ltpf_active = 0;         move16();
    }
    ELSE
    {
        /* Input/Output buffers */
        x = old_x + old_x_len;
        y = old_y + old_y_len;

        /* Input */

        basop_memmove(x, x_in, (L_frame) * sizeof(Word16));

        /* Scaling */
        s0     = sub(s_min(getScaleFactor16_0(old_x, old_x_len), getScaleFactor16_0(old_y, old_y_len)), 1);
        *old_e = sub(*old_e, s0); move16();
        s1     = sub(getScaleFactor16(x, L_frame), 1);
        *x_e   = sub(*x_e, s1); move16();
        s      = sub(*old_e, *x_e);
        IF (s > 0)
        {
            Scale_sig(x, L_frame, sub(s1, s));
            Scale_sig(old_x, old_x_len, s0);
            Scale_sig(old_y, old_y_len, s0);
            *x_e = *old_e; move16();
        }
        ELSE
        {
            Scale_sig(x, L_frame, s1);
            Scale_sig(old_x, old_x_len, add(s0, s));
            Scale_sig(old_y, old_y_len, add(s0, s));
            *old_e = *x_e; move16();
        }

        /* Filtering */
        IF (ltpf_active == 0)
        {
            ltpf_synth_filter(y, x, L_frame / 4, *old_pitch_int, *old_pitch_fr, *old_gain, *old_scale_fac_idx, fs_idx,
                              -1);
        }
        ELSE IF (*mem_ltpf_active == 0)
        {
            ltpf_synth_filter(y, x, L_frame / 4, pitch_int, pitch_fr, gain, scale_fac_idx, fs_idx, 1);
        }
        ELSE IF (sub(pitch_int, *old_pitch_int) == 0 && sub(*old_pitch_fr, pitch_fr) == 0)
        {
            ltpf_synth_filter(y, x, L_frame / 4, pitch_int, pitch_fr, gain, scale_fac_idx, fs_idx, 0);
        }
        ELSE
        {
            ltpf_synth_filter(y, x, L_frame / 4, *old_pitch_int, *old_pitch_fr, *old_gain, *old_scale_fac_idx, fs_idx,
                              -1);
            basop_memmove(z, y - tilt_filter_len[fs_idx], (L_frame / 4 + tilt_filter_len[fs_idx]) * sizeof(Word16));
            ltpf_synth_filter(y, z + tilt_filter_len[fs_idx], L_frame / 4, pitch_int, pitch_fr, gain, scale_fac_idx,
                              fs_idx, 1);
        }
        IF (ltpf_active > 0)
        {
            ltpf_synth_filter(y + L_frame / 4, x + L_frame / 4, 3 * L_frame / 4, pitch_int, pitch_fr, gain,
                              scale_fac_idx, fs_idx, 0);
        }
        ELSE
        {
            basop_memmove(&y[L_frame / 4], &x[L_frame / 4], (3 * L_frame / 4) * sizeof(Word16));
        }

        /* Output */

        basop_memmove(y_out, y, (L_frame) * sizeof(Word16));

        /* Update */

        basop_memmove(old_x, &old_x[L_frame], (old_x_len) * sizeof(Word16));

        basop_memmove(old_y, &old_y[L_frame], (old_y_len) * sizeof(Word16));

        *old_pitch_int   = pitch_int;   move16();
        *old_pitch_fr    = pitch_fr;    move16();
        *old_gain        = gain;        move16();
        *mem_ltpf_active = ltpf_active; move16();
    }

    *old_scale_fac_idx = scale_fac_idx; move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


static void ltpf_synth_filter(Word16 *synth_ltp, Word16 *synth, Word16 length, Word16 pitch_int, Word16 pitch_fr,
                              Word16 gain, Word16 scale_fac_idx, Word16 fs_idx,
                              Word16 fade /* 0=normal, +1=fade-in, -1=fade-out */)
{
    Word16 *x0;
    Word16 *y0;
    Word32  s;
    Word16  alpha, step;
    Word16  i, k;
    Counter j, l;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("ltpf_synth_filter", sizeof(struct {
                   Word16 *x0;
                   Word16 *y0;
                   Word32  s;
                   Word16  alpha, step;
                   Word16  i, k;
                   Counter j, l;
               }));
#endif

    ASSERT(scale_fac_idx >= 0);

    step  = 0; /* initialize just to avoid compiler warning */
    alpha = 0; /* initialize just to avoid compiler warning */
    x0    = &synth_ltp[-pitch_int + inter_filter_shift[fs_idx]];
    y0    = synth;

    alpha = 0; move16();
    IF (fade != 0)
    {
        if (fade < 0)
        {
            alpha = 0x7FFF; move16();
        }

/* step = 1.f/(float)(length); */
        if (sub(length, 5) == 0)
        {
            step = 6553 /*1.f/5.f Q15*/; move16();
        }
        if (sub(length, 10) == 0)
        {
            step = 3276 /*1.f/10.f Q15*/; move16();
        }
        if (sub(length, 15) == 0)
        {
            step = 2184 /*1.f/15.f Q15*/; move16();
        }
        if (sub(length, 20) == 0)
        {
            step = 1638 /*1.f/20.f Q15*/; move16();
        }
        if (sub(length, 30) == 0)
        {
            step = 1092 /*1.f/30.f Q15*/; move16();
        }
        if (sub(length, 40) == 0)
        {
            step = 819 /*1.f/40.f Q15*/; move16();
        }
        if (sub(length, 60) == 0)
        {
            step = 546 /*1.f/60.f Q15*/; move16();
        }
        if (sub(length, 80) == 0)
        {
            step = 409 /*1.f/80.f Q15*/; move16();
        }
        if (sub(length, 120) == 0)
        {
            step = 273 /*1.f/120.f Q15*/; move16();
        }

        if (fade < 0)
            step = negate(step);
    }

    FOR (j = 0; j < length; j++)
    {
        s = L_mult(x0[0], inter_filter[fs_idx][pitch_fr][0]);
        FOR (l = 1; l < inter_filter_len[fs_idx]; l++)
        {
            s = L_mac(s, x0[-l], inter_filter[fs_idx][pitch_fr][l]);
        }
        FOR (l = 0; l < tilt_filter_len[fs_idx]; l++)
        {
            s = L_msu(s, y0[-l], tilt_filter[fs_idx][scale_fac_idx][l]);
        }

        i = msu_r(s, y0[-l], tilt_filter[fs_idx][scale_fac_idx][l]);

        k = mult_r(gain, i);

        if (fade != 0)
            k = mult_r(k, alpha);

        synth_ltp[j] = add(synth[j], k); move16();

        if (fade != 0)
            alpha = add(alpha, step);

        x0++;
        y0++;
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

