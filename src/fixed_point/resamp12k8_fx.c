/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void process_resamp12k8_fx(Word16 x[], Word16 x_len, Word16 mem_in[], Word16 mem_in_len, Word32 mem_50[],
                           Word16 mem_out[], Word16 mem_out_len, Word16 y[], Word16 *y_len, Word16 fs_idx,
                           Word16 frame_dms, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word16 *      buf;
        Word16        index_int, index_frac, len_12k8;
        Word16        resamp_upfac, resamp_off_int, resamp_off_frac, resamp_delay;
        const Word16 *resamp_filt;
        const Word16 *filt_coeff;
        Word16 *      filt_input;
        Word32        filt_output, mem_50_0, mem_50_1;
        Counter       n, m;
        Word32        L_tmp;
    );

    buf = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * (MAX_LEN + MAX_LEN / 8) bytes */

    resamp_upfac    = resamp_params[fs_idx][0]; move16();
    resamp_delay    = resamp_params[fs_idx][1]; move16();
    resamp_off_int  = resamp_params[fs_idx][2]; move16();
    resamp_off_frac = resamp_params[fs_idx][3]; move16();
    resamp_filt     = resamp_filts[fs_idx];     move16();

    len_12k8 = LEN_12K8 / 4 * (frame_dms / 25); move16();
    *y_len   = len_12k8;                        move16();

    /* Init Input Buffer */
    basop_memmove(buf, mem_in, mem_in_len * sizeof(Word16));
    basop_memmove(&buf[mem_in_len], x, x_len * sizeof(Word16));
    basop_memmove(mem_in, &buf[x_len], mem_in_len * sizeof(Word16));

    /* Init Input Indices */
    index_int  = 1; move16();
    index_frac = 0; move16();

    /* Resampling */
    FOR (n = 0; n < len_12k8; n++)
    {
        /* Init Filtering */
        filt_input = &buf[index_int];
        filt_coeff = &resamp_filt[index_frac * resamp_delay * 2];

/* Perform Filtering */
        filt_output = L_mult0(*filt_input, *filt_coeff);
        FOR (m = 1; m < resamp_delay * 2; m++)
        {
            filt_coeff++;
            filt_input++;
            if (*filt_coeff)
            {
                filt_output = L_mac0(filt_output, *filt_input, *filt_coeff);
            }
        }
        y[n] = round_fx(filt_output); move16();

        /* Update Input Indices */
        index_int  = add(index_int, resamp_off_int);
        index_frac = add(index_frac, resamp_off_frac);
        IF (sub(resamp_upfac, index_frac) <= 0)
        {
            index_int  = add(index_int, 1);
            index_frac = sub(index_frac, resamp_upfac);
        }
    }

    /* High Pass Filtering (-3dB at 50Hz) */
    mem_50_0 = mem_50[0]; move32();
    mem_50_1 = mem_50[1]; move32();

    FOR (n = 0; n < len_12k8; n++)
    {
        filt_output = L_mac0(mem_50_0, highpass50_filt_num[0], y[n]);
        L_tmp       = L_mac0(Mpy_32_16(filt_output, highpass50_filt_den[0]), highpass50_filt_num[1], y[n]);
        mem_50_0    = L_add(mem_50_1, L_shl_pos(L_tmp, 1));
        mem_50_1    = L_mac0(Mpy_32_16(filt_output, highpass50_filt_den[1]), highpass50_filt_num[2], y[n]);
        y[n]        = round_fx(filt_output); move16();
    }
    mem_50[0] = mem_50_0; move32();
    mem_50[1] = mem_50_1; move32();

/* Output Buffer */
    basop_memmove(buf, mem_out, mem_out_len * sizeof(Word16));
    basop_memmove(&buf[mem_out_len], y, len_12k8 * sizeof(Word16));
    basop_memmove(y, buf, (*y_len + 1) * sizeof(Word16));
    basop_memmove(mem_out, &buf[len_12k8], mem_out_len * sizeof(Word16));

    Dyn_Mem_Deluxe_Out();
}

