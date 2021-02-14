/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processEncoderEntropy_fl(LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT numbytes, LC3_INT bw_cutoff_bits,
                              LC3_INT bw_cutoff_idx, LC3_INT lastnz, LC3_INT N, LC3_INT lsbMode, LC3_INT gg_idx, LC3_INT num_tns_filters,
                              LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* scf_idx, LC3_INT fac_ns_idx)
{
    LC3_UINT8* ptr;
    LC3_INT      i = 0, submodeMSB = 0, submodeLSB = 0, tmp = 0, gainMSB = 0, gainLSB = 0;
    LC3_INT      gainMSBbits[4] = {1, 1, 2, 2}, gainLSBbits[4] = {0, 1, 0, 1};

    *bp_side   = numbytes - 1;
    *mask_side = 1;
    ptr        = bytes;

    /* Bandwidth */
    if (bw_cutoff_bits > 0) {
        write_uint_backward_fl(ptr, bp_side, mask_side, bw_cutoff_idx, bw_cutoff_bits);
    }

    /* Last non zero touple */
    write_uint_backward_fl(ptr, bp_side, mask_side, lastnz / 2 - 1, ceil(LC3_LOG2(N / 2)));

    /* LSB mode bit */
    write_bit_backward_fl(ptr, bp_side, mask_side, lsbMode);

    /* Global gain */
    write_uint_backward_fl(ptr, bp_side, mask_side, gg_idx, 8);

    /* TNS activation flag */
    for (i = 0; i < num_tns_filters; i++) {
        write_bit_backward_fl(ptr, bp_side, mask_side, MIN(1, tns_order[i]));
    }

    /* LTPF activation flag */
    write_bit_backward_fl(ptr, bp_side, mask_side, ltpf_idx[0]);

    /* SNS-VQ 1st stage */
    write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[0], 5);
    write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[1], 5);

    /* SNS-VQ 2nd stage side-info (3-4 bits) */
    submodeMSB = scf_idx[2] / 2;
    submodeLSB = scf_idx[2] & 1;
    write_bit_backward_fl(ptr, bp_side, mask_side, submodeMSB);
    gainMSB = scf_idx[3] >> (gainLSBbits[scf_idx[2]]);
    gainLSB = scf_idx[3] & 1;
    write_uint_backward_fl(ptr, bp_side, mask_side, gainMSB, gainMSBbits[scf_idx[2]]);
    write_bit_backward_fl(ptr, bp_side, mask_side, scf_idx[4]);

    /* SNS-VQ 2nd stage MPVQ data (24-25 bits) */
    if (submodeMSB == 0) {
        if (submodeLSB == 0) {
            tmp = scf_idx[6] + 2;
        } else {
            tmp = gainLSB;
        }

        tmp = tmp * 2390004 + scf_idx[5];
        write_uint_backward_fl(ptr, bp_side, mask_side, tmp, 25);
    } else {
        tmp = scf_idx[5];

        if (submodeLSB != 0) {
            tmp = 2 * tmp + gainLSB + 15158272;
        }

        write_uint_backward_fl(ptr, bp_side, mask_side, tmp, 24);
    }

    /* LTPF data */
    if (ltpf_idx[0] == 1) {
        write_uint_backward_fl(ptr, bp_side, mask_side, ltpf_idx[1], 1);
        write_uint_backward_fl(ptr, bp_side, mask_side, ltpf_idx[2], 9);
    }

    /* Noise factor */
    write_uint_backward_fl(ptr, bp_side, mask_side, fac_ns_idx, 3);
}

void write_uint_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT val, LC3_INT numbits)
{
    LC3_INT k = 0, bit = 0;

    for (k = 0; k < numbits; k++) {
        bit = val & 1;
        write_bit_backward_fl(ptr, bp_side, mask_side, bit);
        val = val / 2;
    }
}

void write_bit_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT bit)
{
    if (bit == 0) {
        ptr[*bp_side] = ptr[*bp_side] & (255 - *mask_side);
    } else {
        ptr[*bp_side] = ptr[*bp_side] | *mask_side;
    }

    if (*mask_side == 128) {
        *mask_side = 1;
        *bp_side   = *bp_side - 1;
    } else {
        *mask_side = *mask_side * 2;
    }
}
