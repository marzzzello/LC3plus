/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


static void   Parcor2Index(const Word16 parCoeff[] /*Q15*/, Word16 index[], Word16 order);
static void   Index2Parcor(const Word16 index[], Word16 parCoeff[], Word16 order);
static Word32 FIRLattice(Word16 order, const Word16 *parCoeff /*Q15*/, Word32 *state, Word32 x /* Q0 */);

/*************************************************************************/

void processTnsCoder_fx(Word16 *bits, Word16 indexes[], Word32 x[], Word16 BW_cutoff_idx, Word16 order[],
                        Word16 *numfilters, Word16 enable_lpc_weighting, Word16 nSubdivisions, Word16 frame_dms,
                        Word16 maxLen, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word16 *      tmpbuf;
        Word32 *      rxx, epsP, *state, L_tmp, *A, predictionGain, alpha;
        Word16 *      RC, inv;
        Word16        n, n2, headroom, shift, tmp, shifts, facs, facs_e, stopfreq, xLen, maxOrder;
        Word16        startfreq[TNS_NUMFILTERS_MAX];
        const Word16 *subdiv_startfreq, *subdiv_stopfreq;
        Counter       i, j, iSubdivisions, lag;
        Word8 *       LevinsonBuffer;
    );

    /* Buffer alignment */
    tmpbuf = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN */

    rxx = (Word32 *)scratchAlign(tmpbuf, sizeof(*tmpbuf) * maxLen); /* Size = 4 * (MAXLAG + 1) = 36 bytes */

    state = (Word32 *)scratchAlign(rxx, sizeof(*rxx) * (MAXLAG + 1)); /* Size = 4 * MAXLAG = 32 bytes */

    A = (Word32 *)scratchAlign(state, sizeof(*state) * MAXLAG); /* Size = 4 * (MAXLAG + 1) = 36 bytes */

    RC = (Word16 *)scratchAlign(A, sizeof(*A) * (MAXLAG + 1)); /* Size = 2 * MAXLAG = 16 bytes */

    LevinsonBuffer = (Word8 *)scratchAlign(RC, sizeof(*RC) * (MAXLAG)); /* Size = 4 * (M_LTPF + 1) = 100 bytes */

    /* Init */
    *bits       = 0;                                move16();
    maxOrder    = MAXLAG;                           move16();
    *numfilters = 1;                                move16();
    xLen        = BW_cutoff_bin_all[BW_cutoff_idx]; move16();

    SWITCH (frame_dms)
    {
    case 25:
        startfreq[0]     = 3;                                         move16();
        subdiv_startfreq = tns_subdiv_startfreq_2_5ms[BW_cutoff_idx]; move16();
        subdiv_stopfreq  = tns_subdiv_stopfreq_2_5ms[BW_cutoff_idx];  move16();
        xLen             = shr_pos(xLen, 2);
        maxOrder         = 4; move16();
        BREAK;
    case 50:
        startfreq[0]     = 6;                                       move16();
        subdiv_startfreq = tns_subdiv_startfreq_5ms[BW_cutoff_idx]; move16();
        subdiv_stopfreq  = tns_subdiv_stopfreq_5ms[BW_cutoff_idx];  move16();
        xLen             = shr_pos(xLen, 1);
        maxOrder         = 4;
        BREAK;
    default:                                                    /* 100 */
        startfreq[0]     = 12;                                  move16();
        subdiv_startfreq = tns_subdiv_startfreq[BW_cutoff_idx]; move16();
        subdiv_stopfreq  = tns_subdiv_stopfreq[BW_cutoff_idx];  move16();
        BREAK;
    }

    IF (sub(BW_cutoff_idx, 3) >= 0 && frame_dms >= 50)
    {
        *numfilters  = 2;
        startfreq[1] = shr_pos(xLen, 1);
    }

    basop_memset(state, 0, MAXLAG * sizeof(*state));

    FOR (j = 0; j < *numfilters; j++)
    {
        basop_memset(rxx, 0, (maxOrder + 1) * sizeof(*rxx));

        FOR (iSubdivisions = 0; iSubdivisions < nSubdivisions; iSubdivisions++)
        {
            n = sub(subdiv_stopfreq[nSubdivisions * j + iSubdivisions],
                    subdiv_startfreq[nSubdivisions * j + iSubdivisions]);

            /*norms[iFilter][iSubdivisions] = norm2FLOAT(pSpectrum+iStartLine, iEndLine-iStartLine);*/
            headroom = getScaleFactor32(x + subdiv_startfreq[nSubdivisions * j + iSubdivisions], n);

            /* Calculate norm of spectrum band */
            L_tmp = Norm32Norm(x + subdiv_startfreq[nSubdivisions * j + iSubdivisions], headroom, n, &shift);

            /* Rounding to avoid overflow when computing the autocorrelation below */
            tmp   = sub(norm_l(L_tmp), 1);
            L_tmp = L_shl(L_tmp, tmp);
            shift = sub(shift, tmp);
            L_tmp = L_add(L_tmp, 0x8000);
            L_tmp = L_and(L_tmp, 0x7FFF0000);

            IF (L_tmp == 0)
            {
                rxx[0] = 0x7FFFFFFF; move32();
                basop_memset(&rxx[1], 0, (maxOrder) * sizeof(*rxx));
                BREAK;
            }

            /* get pre-shift for autocorrelation */
            tmp    = sub(shift, norm_l(L_tmp)); /* exponent for normalized L_tmp */
            tmp    = shr_pos(sub(1, tmp), 1);   /* pre-shift to apply before autocorrelation */
            shifts = s_min(tmp, headroom);

            /* calc normalization factor */
            facs_e = shl_pos(sub(tmp, shifts), 1);

            SWITCH (frame_dms)
            {
            case 25: facs_e = add(facs_e, 1); BREAK;
            case 50: facs_e = add(facs_e, 1); BREAK;
            case 100: BREAK;
            }

            tmp   = sub(1, shl_pos(tmp, 1));       /* exponent of autocorrelation */
            L_tmp = L_shl(L_tmp, sub(shift, tmp)); /* shift L_tmp to that exponent */
            /* calc factor (with 2 bits headroom for sum of 3 subdivisions) */
            facs = div_s(0x2000, round_fx(L_tmp)); /* L_tmp is >= 0x2000000 */

            FOR (i = 0; i < n; i++)
            {
                tmpbuf[i] = round_fx_sat(
                    L_shl_sat(x[subdiv_startfreq[nSubdivisions * j + iSubdivisions] + i], shifts)); move16();
            }

            FOR (lag = 0; lag <= maxOrder; lag++)
            {
                n2 = sub(n, lag);
                L_tmp = L_deposit_l(0);
                FOR (i = 0; i < n2; i++)
                {
                    L_tmp = L_mac0(L_tmp, tmpbuf[i], tmpbuf[i + lag]);
                }
                if (lag != 0)
                    L_tmp = Mpy_32_32(L_tmp, tnsAcfWindow[lag - 1]);

                L_tmp = Mpy_32_16(L_tmp, facs);
                L_tmp = L_shl(L_tmp, facs_e);

                rxx[lag] = L_add(rxx[lag], L_tmp); move32();
            }
        }

        /* Levinson-Durbin */
        processLevinson_fx(A, rxx, maxOrder, RC, &epsP, LevinsonBuffer);

        /* Prediction Gain */
        shift          = norm_l(epsP);
        inv            = div_s(16383, extract_h(L_shl_pos(epsP, shift)));
        predictionGain = Mpy_32_32(rxx[0], Mpy_32_16(L_sub(MAX_32, Mpy_32_16(L_shl(epsP, shift), inv)), inv));

        IF (L_sub(predictionGain, L_shr_pos_pos(0x30000000, shift)) > 0)
        {
            /* If Prediction Gain is low */
            test();
            IF (enable_lpc_weighting != 0 && L_sub(predictionGain, L_shr_pos_pos(0x40000000, shift)) < 0)
            {
                /* LPC weighting */
                alpha = L_add(0x6CCCCCCD,
                              Mpy_32_32(0x13333333, L_shl_pos(L_sub(L_shl_pos(predictionGain, shift), 0x30000000), 3)));
                L_tmp = alpha;
                FOR (i = 1; i < maxOrder; i++)
                {
                    A[i]  = Mpy_32_32(A[i], L_tmp); move32();
                    L_tmp = Mpy_32_32(L_tmp, alpha);
                }
                A[maxOrder] = Mpy_32_32(A[maxOrder], L_tmp); move32();

                /* LPC -> RC */
                lpc2rc(A, RC, maxOrder);
            }

            /* Reflection Coefficients Quantization */
            Parcor2Index(RC, &indexes[MAXLAG * j], maxOrder);

            /* reduce filter order by truncating trailing zeros */
            i = sub(maxOrder, 1);
            WHILE ((i >= 0) && (indexes[MAXLAG * j + i] == INDEX_SHIFT))
            {
                i = sub(i, 1);
            }
            order[j] = add(i, 1);

            /* Count bits */
            L_tmp = L_deposit_l(ac_tns_order_bits[enable_lpc_weighting][order[j] - 1]);
            FOR (i = 0; i < order[j]; i++)
            {
                L_tmp = L_add(L_tmp, L_deposit_l(ac_tns_coef_bits[i][indexes[MAXLAG * j + i]]));
            }
            *bits = add(*bits, add(2, extract_l(L_shr_pos(L_sub(L_tmp, 1), 11)))); move16();

            /* Unquantize Reflection Coefficients */
            Index2Parcor(&indexes[MAXLAG * j], RC, order[j]);

            /* Stop frequency */
            stopfreq = xLen; move16();
            IF (sub(*numfilters, 2) == 0 && j == 0)
            {
                stopfreq = startfreq[1];
            }

            /* Filter */
            FOR (i = startfreq[j]; i < stopfreq; i++)
            {
                x[i] = FIRLattice(order[j], RC, state, x[i]); move32();
            }
        }
        ELSE
        {
            /* TNS disabled */
            *bits    = add(*bits, 1);
            order[j] = 0;
        }
    }

    Dyn_Mem_Deluxe_Out();
}

/*************************************************************************/

static void Parcor2Index(const Word16 parCoeff[] /*Q15*/, Word16 index[], Word16 order)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word16  iIndex;
        Word16  x;
    );

    FOR (i = 0; i < order; i++)
    {
        move16(); move16();
        iIndex = 1;
        x      = parCoeff[i];

        WHILE ((iIndex < TNS_COEF_RES) && (x > tnsQuantThr[iIndex - 1]))
        {
            iIndex = add(iIndex, 1);
        }
        index[i] = sub(iIndex, 1); move16();
    }

    Dyn_Mem_Deluxe_Out();
}

static void Index2Parcor(const Word16 index[], Word16 parCoeff[], Word16 order)
{
    Counter i;
    FOR (i = 0; i < order; i++)
    {
        parCoeff[i] = tnsQuantPts[index[i]]; move16();
    }
}

static Word32 FIRLattice(Word16 order, const Word16 *parCoeff /*Q15*/, Word32 *state, Word32 x /* Q0 */)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32  tmpSave, tmp;
    );

    tmpSave = L_add(x, 0);

    FOR (i = 0; i < order - 1; i++)
    {
        tmp      = L_add(state[i], Mpy_32_16(x, parCoeff[i]));
        x        = L_add(x, Mpy_32_16(state[i], parCoeff[i])); /* exponent: 31+0 */
        state[i] = tmpSave;                                    move32();
        tmpSave  = L_add(tmp, 0);
    }

    /* last stage: only need half operations */
    x                = L_add(x, Mpy_32_16(state[order - 1], parCoeff[order - 1]));
    state[order - 1] = tmpSave; move32();
    Dyn_Mem_Deluxe_Out();
    return x;
}

