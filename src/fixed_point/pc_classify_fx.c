/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"

#include "constants.h"
#include "functions.h"


#define BLOCK_SIZE 3
#define THR1 8
#define FAC 9830 /* 0.3 */

void peakDetector_fx(Word16 in_sig[], Word16 yLen, Word16 *xover);

void processPCclassify_fx(Word16 pitch_present, Word16 frame_dms, Word16 q_old_d_fx[], Word16 q_old_res_fx[],
                          Word16 yLen, Word16 spec_inv_idx, Word16 stab_fac, Word16 prev_bfi, Word16 *bfi)
{
    Dyn_Mem_Deluxe_In(
        Word16  maxPitchBin, xover;
        Counter i;
        Word16  s, tmp16, full_nrg16, part_nrg16;
        Word32  full_nrg, part_nrg;
    );

    IF (sub(prev_bfi, 1) == 0)
    {
        *bfi = 1;
    }
    /* Apply classifier only if lower than 2khz signal */
    ELSE IF (sub(i_mult(spec_inv_idx, 10), shl_pos(frame_dms, 2)) < 0 )
    {
        IF (sub(stab_fac, 16384 /* 0.5 */) < 0)
        {
            *bfi = 1;
        }
        ELSE IF (sub(pitch_present, 1) == 0)
        {
            maxPitchBin = 8;  move16();
            IF (sub(frame_dms, 50) == 0)
            {
                maxPitchBin = 4;  move16();
            }

            /* avoid phase discontinuity in low frequencies */
            peakDetector_fx(q_old_d_fx, yLen, &xover);
            test();
            IF (sub(spec_inv_idx, xover) < 0 || sub(spec_inv_idx, maxPitchBin) < 0)
            {
                *bfi = 1;
            }
        }
        ELSE
        {
            s = getScaleFactor16(q_old_res_fx, yLen);

            part_nrg = 0;  move32();
            FOR (i = 0; i < spec_inv_idx; i++)
            {
                tmp16    = shl_sat(q_old_res_fx[i], sub(s, 4));
                part_nrg = L_mac0(part_nrg, tmp16, tmp16); /* exp = 2s - 8 */
            }

            full_nrg = part_nrg;  move32();
            FOR (i = spec_inv_idx; i < yLen; i++)
            {
                tmp16    = shl_sat(q_old_res_fx[i], sub(s, 4));
                full_nrg = L_mac0(full_nrg, tmp16, tmp16); /* exp = 2s - 8 */
            }

            s          = getScaleFactor32(&full_nrg, 1);
            full_nrg16 = extract_h(L_shl(full_nrg, s));
            part_nrg16 = extract_h(L_shl(part_nrg, s));

            tmp16 = mult(full_nrg16, 9830 /* 0.3 */);

            IF (part_nrg16 < tmp16)
            {
                *bfi = 1;
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
}

void peakDetector_fx(Word16 in_sig[], Word16 yLen, Word16 *xover)
{
    Dyn_Mem_Deluxe_In(
        Counter i, j;
        Word16  tmp16, c, s, s2, mean_block_nrg16;
        Word32  maxPeak, tmp32;
        Word32  mean_block_nrg, block_cent;
        Word16  cur_max, prev_max, next_max;
    );

    *xover = 0;

    s = getScaleFactor16(in_sig, yLen);

    mean_block_nrg = 0;  move32();
    FOR (i = 0; i < yLen; i++)
    {
        tmp16          = shl_sat(in_sig[i], sub(s, 4));
        mean_block_nrg = L_mac0(mean_block_nrg, tmp16, tmp16); /* exp = 2s - 8 */
    }

    s2               = getScaleFactor16(&yLen, 1);
    c                = shl(yLen, s2);
    mean_block_nrg16 = div_l(mean_block_nrg, c);                                        /* exp = 2s - 8 - s2 - 1 */
    mean_block_nrg   = L_shl(L_mult0(mean_block_nrg16, BLOCK_SIZE * THR1), add(4, s2)); /* exp = 2s - 5 */

    maxPeak = 0;  move32();
    c = sub(yLen, 2 * BLOCK_SIZE);

    test();
    IF (abs_s(in_sig[0]) >= abs_s(in_sig[1]))
    {
        block_cent = 0;  move32();
        FOR (j = 0; j <= 1; j++)
        {
            tmp16      = shl_sat(in_sig[j], sub(s, 2));
            block_cent = L_mac0(block_cent, tmp16, tmp16); /* -> exp = 2s - 4 */
        }
        block_cent = L_shr(block_cent, 1); /* exp = 2s - 5 */

        IF (L_sub(block_cent, mean_block_nrg) > 0)
        {
            cur_max = abs_s(in_sig[0]);
            cur_max = MAX(abs_s(in_sig[1]), cur_max);

            next_max = abs_s(in_sig[-1 + BLOCK_SIZE]);
            next_max = MAX(abs_s(in_sig[-1 + BLOCK_SIZE + 1]), next_max);
            next_max = MAX(abs_s(in_sig[-1 + BLOCK_SIZE + 2]), next_max);

            IF (sub(cur_max, next_max) > 0)
            {
                maxPeak = block_cent;  move32();
                *xover = 1;
            }
        }
    }

    FOR (i = 0; i < BLOCK_SIZE; i++)
    {
        test();
        IF (abs_s(in_sig[i + 1]) >= abs_s(in_sig[i]) && abs_s(in_sig[i + 1]) >= abs_s(in_sig[i + 2]))
        {
            block_cent = 0;  move32();
            FOR (j = 0; j < BLOCK_SIZE; j++)
            {
                tmp16      = shl_sat(in_sig[i + j], sub(s, 2));
                block_cent = L_mac0(block_cent, tmp16, tmp16); /* -> exp = 2s - 4 */
            }
            block_cent = L_shr(block_cent, 1); /* exp = 2s - 5 */

            IF (L_sub(block_cent, mean_block_nrg) > 0)
            {
                cur_max = abs_s(in_sig[i]);
                cur_max = MAX(abs_s(in_sig[i + 1]), cur_max);
                cur_max = MAX(abs_s(in_sig[i + 2]), cur_max);

                prev_max = 0;  move16();
                FOR (j = i - BLOCK_SIZE; j <= i - 1; j++)
                {
                    IF (j > 0)
                    {
                        prev_max = MAX(abs_s(in_sig[j]), prev_max);
                    }
                }

                next_max = abs_s(in_sig[i + BLOCK_SIZE]);
                next_max = MAX(abs_s(in_sig[i + BLOCK_SIZE + 1]), next_max);
                next_max = MAX(abs_s(in_sig[i + BLOCK_SIZE + 2]), next_max);

                test();
                IF (sub(cur_max, prev_max) >= 0 && sub(cur_max, next_max) > 0)
                {
                    IF (L_sub(block_cent, maxPeak) >= 0)
                    {
                        maxPeak = block_cent;  move32();
                        *xover = sub(add(i, BLOCK_SIZE), 1);
                    }
                    ELSE
                    {
                        tmp32 = L_mult(FAC, extract_h(maxPeak));

                        tmp16 = extract_l(L_shr(maxPeak, 1));
                        tmp16 = s_and(tmp16, 0x7fff);
                        tmp16 = mult(FAC, tmp16);
                        tmp32 = L_add_sat(tmp32, tmp16);

                        IF (L_sub(block_cent, tmp32) > 0)
                        {
                            *xover = sub(add(i, BLOCK_SIZE), 1);
                        }
                    }
                }
            }
        }
    }

    FOR (i = BLOCK_SIZE; i <= c; i++)
    {
        test();
        IF (abs_s(in_sig[i + 1]) >= abs_s(in_sig[i]) && abs_s(in_sig[i + 1]) >= abs_s(in_sig[i + 2]))
        {
            block_cent = 0;  move32();
            FOR (j = 0; j < BLOCK_SIZE; j++)
            {
                tmp16      = shl_sat(in_sig[i + j], sub(s, 2));
                block_cent = L_mac0(block_cent, tmp16, tmp16); /* -> exp = 2s - 4 */
            }
            block_cent = L_shr(block_cent, 1); /* exp = 2s - 5 */

            IF (L_sub(block_cent, mean_block_nrg) > 0)
            {
                cur_max = abs_s(in_sig[i]);
                cur_max = MAX(abs_s(in_sig[i + 1]), cur_max);
                cur_max = MAX(abs_s(in_sig[i + 2]), cur_max);

                prev_max = abs_s(in_sig[i - BLOCK_SIZE]);
                prev_max = MAX(abs_s(in_sig[i - BLOCK_SIZE + 1]), prev_max);
                prev_max = MAX(abs_s(in_sig[i - BLOCK_SIZE + 2]), prev_max);

                next_max = abs_s(in_sig[i + BLOCK_SIZE]);
                next_max = MAX(abs_s(in_sig[i + BLOCK_SIZE + 1]), next_max);
                next_max = MAX(abs_s(in_sig[i + BLOCK_SIZE + 2]), next_max);

                test();
                IF (sub(cur_max, prev_max) >= 0 && sub(cur_max, next_max) > 0)
                {
                    IF (L_sub(block_cent, maxPeak) >= 0)
                    {
                        maxPeak = block_cent;  move32();
                        *xover = sub(add(i, BLOCK_SIZE), 1);
                    }
                    ELSE
                    {
                        tmp32 = L_mult(FAC, extract_h(maxPeak));

                        tmp16 = extract_l(L_shr(maxPeak, 1));
                        tmp16 = s_and(tmp16, 0x7fff);
                        tmp16 = mult(FAC, tmp16);
                        tmp32 = L_add_sat(tmp32, tmp16);

                        IF (L_sub(block_cent, tmp32) > 0)
                        {
                            *xover = sub(add(i, BLOCK_SIZE), 1);
                        }
                    }
                }
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
}


