/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit);
static void read_uint_fl(LC3_INT nbits, LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* val);

void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit)
{
    if (ptr[*bp_side] & *mask_side) {
        *bit = 1;
    } else {
        *bit = 0;
    }

    if (*mask_side == 128) {
        *mask_side = 1;
        *bp_side   = *bp_side - 1;
    } else {
        *mask_side = *mask_side * 2;
    }
}

void read_uint_fl(LC3_INT nbits, LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* val)
{
    LC3_INT bit = 0, i = 0;

    read_bit_fl(ptr, mask_side, bp_side, val);

    for (i = 1; i < nbits; i++) {
        read_bit_fl(ptr, mask_side, bp_side, &bit);
        *val = *val + (bit << i);
    }
}

#ifdef ENABLE_PADDING
LC3_INT paddingDec_fl(LC3_UINT8* bytes, LC3_INT numbytes, LC3_INT N, LC3_INT bw_cutoff_bits, LC3_INT* total_padding)
{
    LC3_INT lastnz_threshold;
    LC3_INT val, padding_len_bits, padding_len;

    LC3_INT bp_side;

    LC3_INT    mask_side;
    LC3_UINT8* ptr = bytes;

    LC3_INT lastnz;
    LC3_INT bw_cutoff_idx;
    LC3_INT nbits = ceil(LC3_LOG2(N / 2));

    *total_padding = 0;

    bp_side   = numbytes - 1; /* Matlab offset by 1 */
    mask_side = 1;

    if (bp_side < 19 || bp_side >= LC3_MAX_BYTES) {
        return 1;
    }
    
    ptr = bytes;

    if (bw_cutoff_bits > 0) {
        read_uint_fl(bw_cutoff_bits, ptr, &mask_side, &bp_side, &bw_cutoff_idx);
    }

    read_uint_fl(nbits, ptr, &mask_side, &bp_side, &lastnz);

    lastnz_threshold = (1 << nbits) - 1 - 1;

    while (lastnz == lastnz_threshold) {
        padding_len_bits = 16 - nbits - bw_cutoff_bits - 4;

        /*Read padding length*/
        read_uint_fl(padding_len_bits, ptr, &mask_side, &bp_side, &padding_len);

        /* Read 4 reserved bits */
        read_uint_fl(4, ptr, &mask_side, &bp_side, &val);

        /* Discard padding length bytes */
        bp_side = bp_side - padding_len;

        *total_padding = *total_padding + padding_len + 2;

        /* check if minimum payload size is reached */
        if ((numbytes - *total_padding) < 20) {
            return 1;
        }

        /* Read bandwidth bits */
        if (bw_cutoff_bits > 0) {
            read_uint_fl(bw_cutoff_bits, ptr, &mask_side, &bp_side, &bw_cutoff_idx);
        }

        read_uint_fl(nbits, ptr, &mask_side, &bp_side, &lastnz);
    }
    return 0;
}
#endif

void processDecoderEntropy_fl(LC3_UINT8* bytes, LC3_INT numbytes, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT N, LC3_INT fs_idx,
                              LC3_INT bw_cutoff_bits, LC3_INT* bfi, LC3_INT* gg_idx, LC3_INT* scf_idx, LC3_INT* fac_ns_idx,
                              LC3_INT* tns_numfilters, LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* bw_cutoff_idx, LC3_INT* lastnz,
                              LC3_INT* lsbMode, LC3_INT frame_dms)
{

    LC3_INT plc_trigger_bw = 0, plc_trigger_last_nz = 0, plc_trigger_SNS1 = 0, plc_trigger_SNS2 = 0, tmp = 0, bit = 0,
        submodeMSB = 0, i = 0, ltpf_tmp[3] = {0}, ind = 0, submodeLSB = 0;
    LC3_UINT8* ptr;
    LC3_INT      gainMSBbits[4] = {1, 1, 2, 2};

    *bp_side   = numbytes - 1; /* Matlab offset by 1 */
    *mask_side = 1;
    ptr        = bytes;

    plc_trigger_bw      = 1; /* Bandwidth */
    plc_trigger_last_nz = 1; /* Last non-zero tuple */
    plc_trigger_SNS1    = 1; /* SNS-VQ 2nd stage MPVQ data (24-25 bits) */
    plc_trigger_SNS2    = 1; /* SNS-VQ 2nd stage MPVQ data (24-25 bits) */

    /* Bandwidth */
    if (bw_cutoff_bits > 0) {
        read_uint_fl(bw_cutoff_bits, ptr, mask_side, bp_side, bw_cutoff_idx);

        if (fs_idx < *bw_cutoff_idx) {
            *bfi = plc_trigger_bw;

            if (*bfi) {
                return;
            }
        }
    } else {
    	*bw_cutoff_idx = fs_idx;
    }

    /* Number of TNS filters */
    if (*bw_cutoff_idx < 3 || frame_dms == 25) {
        *tns_numfilters = 1;
    } else {
        *tns_numfilters = 2;
    }

    /* Last non-zero tuple */
    read_uint_fl(ceil(LC3_LOG2(N / 2)), ptr, mask_side, bp_side, lastnz);
    *lastnz = (*lastnz + 1) * 2;

    if (*lastnz > N) {
        *bfi = plc_trigger_last_nz;
        if (*bfi) {
            return;
        }
    }

    /* LSB mode bit */
    read_bit_fl(ptr, mask_side, bp_side, lsbMode);

    /* Global gain */
    read_uint_fl(8, ptr, mask_side, bp_side, gg_idx);

    /* TNS activation flag */
    for (i = 0; i < *tns_numfilters; i++) {
        read_bit_fl(ptr, mask_side, bp_side, &bit);
        tns_order[i] = bit;
    }

    /* LTPF activation flag */
    read_bit_fl(ptr, mask_side, bp_side, &ltpf_tmp[0]);

    /* SNS-VQ 1st stage */
    read_uint_fl(5, ptr, mask_side, bp_side, &scf_idx[0]);
    read_uint_fl(5, ptr, mask_side, bp_side, &scf_idx[1]);

    /* SNS-VQ 2nd stage side-info (3-4 bits) */
    read_bit_fl(ptr, mask_side, bp_side, &submodeMSB);
    scf_idx[2] = submodeMSB * 2;

    read_uint_fl(gainMSBbits[scf_idx[2]], ptr, mask_side, bp_side, &scf_idx[3]);
    read_bit_fl(ptr, mask_side, bp_side, &scf_idx[4]);

    /* SNS-VQ 2nd stage MPVQ data (24-25 bits) */
    if (submodeMSB == 0) {
        read_uint_fl(25, ptr, mask_side, bp_side, &tmp);
        if (tmp >= 33460056) {
            *bfi = plc_trigger_SNS1;
            if (*bfi) {
                return;
            }
        }

        ind        = floor(tmp / 2390004);
        scf_idx[5] = tmp - ind * 2390004;

        if (ind < 2) {
            submodeLSB = 1;
            scf_idx[3] = scf_idx[3] * 2 + ind;
            scf_idx[6] = -2;
        } else {
            submodeLSB = 0;
            scf_idx[6] = ind - 2;
        }

    } else {
        read_uint_fl(24, ptr, mask_side, bp_side, &tmp);

        if (tmp >= 16708096) {
            *bfi = plc_trigger_SNS2;
            if (*bfi) {
                return;
            }
        }

        if (tmp >= 15158272) {
            submodeLSB = 1;
            scf_idx[3] = scf_idx[3] * 2 + ((tmp - 15158272) & 1);
            scf_idx[5] = floor((tmp - 15158272) / 2);
            scf_idx[6] = -2;
        } else {
            submodeLSB = 0;
            scf_idx[5] = tmp;
            scf_idx[6] = -1;
        }
    }

    scf_idx[2] = scf_idx[2] + submodeLSB;

    /* LTPF data */
    if (ltpf_tmp[0] == 1) {
        read_bit_fl(ptr, mask_side, bp_side, &ltpf_tmp[1]);
        read_uint_fl(9, ptr, mask_side, bp_side, &ltpf_tmp[2]);
    } else {
        ltpf_tmp[1] = 0;
        ltpf_tmp[2] = 0;
    }

    for (i = 0; i < 3; i++) {
        ltpf_idx[i] = ltpf_tmp[i];
    }

    /* Noise factor */
    read_uint_fl(3, ptr, mask_side, bp_side, fac_ns_idx);
}
