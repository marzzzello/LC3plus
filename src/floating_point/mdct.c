/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef ENABLE_HR_MODE
static const LC3_FLOAT* mdct_window(LC3_INT length, LC3_INT frame_dms, LC3_INT hrmode)
{
    if (frame_dms == 100) {
        switch (length) {
        case 80:
            return MDCT_WINS_10ms[hrmode][0];
        case 160:
            return MDCT_WINS_10ms[hrmode][1];
        case 240:
            return MDCT_WINS_10ms[hrmode][2];
        case 320:
            return MDCT_WINS_10ms[hrmode][3];
        case 480:
            return MDCT_WINS_10ms[hrmode][4];
        case 960:
        	return MDCT_WINS_10ms[hrmode][5];
        default:
            return NULL;
        }
    }
    else if (frame_dms == 50) {
        switch (length) {
        case 40:
            return MDCT_WINS_5ms[hrmode][0];
        case 80:
            return MDCT_WINS_5ms[hrmode][1];
        case 120:
            return MDCT_WINS_5ms[hrmode][2];
        case 160:
            return MDCT_WINS_5ms[hrmode][3];
        case 240:
            return MDCT_WINS_5ms[hrmode][4];
        case 480:
            return MDCT_WINS_5ms[hrmode][5];
        default:
            return NULL;
        }
    }
    else if (frame_dms == 25) {
        switch (length) {
        case 20:
            return MDCT_WINS_2_5ms[hrmode][0];
        case 40:
            return MDCT_WINS_2_5ms[hrmode][1];
        case 60:
            return MDCT_WINS_2_5ms[hrmode][2];
        case 80:
            return MDCT_WINS_2_5ms[hrmode][3];
        case 120:
            return MDCT_WINS_2_5ms[hrmode][4];
        case 240:
        	return MDCT_WINS_2_5ms[hrmode][5];
        default:
            return NULL;
        }
    }
    return NULL;
}
#else
static const LC3_FLOAT* mdct_window(LC3_INT length, LC3_INT frame_dms)
{
    if (frame_dms == 100) {
        switch (length) {
        case 80:
            return MDCT_WINS_10ms[0];
        case 160:
            return MDCT_WINS_10ms[1];
        case 240:
            return MDCT_WINS_10ms[2];
        case 320:
            return MDCT_WINS_10ms[3];
        case 480:
            return MDCT_WINS_10ms[4];
        default:
            return NULL;
        }
    }
    if (frame_dms == 50) {
        switch (length) {
        case 40:
            return MDCT_WINS_5ms[0];
        case 80:
            return MDCT_WINS_5ms[1];
        case 120:
            return MDCT_WINS_5ms[2];
        case 160:
            return MDCT_WINS_5ms[3];
        case 240:
            return MDCT_WINS_5ms[4];
        default:
            return NULL;
        }
    }
    if (frame_dms == 25) {
        switch (length) {
        case 20:
            return MDCT_WINS_2_5ms[0];
        case 40:
            return MDCT_WINS_2_5ms[1];
        case 60:
            return MDCT_WINS_2_5ms[2];
        case 80:
            return MDCT_WINS_2_5ms[3];
        case 120:
            return MDCT_WINS_2_5ms[4];
        default:
            return NULL;
        }
    }
    return NULL;
}
#endif

#ifdef ENABLE_HR_MODE
void mdct_init(Mdct* mdct, LC3_INT length, LC3_INT frame_dms, LC3_INT hrmode)
#else
void mdct_init(Mdct* mdct, LC3_INT length, LC3_INT frame_dms)
#endif
{
    if (frame_dms == 100) {
        mdct->leading_zeros = 3 * length / 8;
    } 
    else if (frame_dms == 50) {
        mdct->leading_zeros = length / 4;
    } 
    else if (frame_dms == 25) {
        mdct->leading_zeros = 0;
    }
    else {
        assert(!"invalid frame_ms");
    }

    mdct->length     = length;
    mdct->mem_length = length - mdct->leading_zeros;
#ifdef ENABLE_HR_MODE
    mdct->window     = mdct_window(length, frame_dms, hrmode);
#else
    mdct->window     = mdct_window(length, frame_dms);
#endif
    mdct->mem        = calloc(sizeof(*mdct->mem), mdct->mem_length);
    dct4_init(&mdct->dct, length);
}

void mdct_free(Mdct* mdct)
{
    if (mdct) {
        free(mdct->mem);
        dct4_free(&mdct->dct);
        memset(mdct, 0, sizeof(*mdct));
    }
}

void mdct_apply(const LC3_FLOAT* input, LC3_FLOAT* output, Mdct* mdct)
{
    LC3_FLOAT tmp[MAX_LEN * 2] = {0};
    LC3_INT   i = 0;

    move_float(tmp, mdct->mem, mdct->mem_length);
    move_float(tmp + mdct->mem_length, input, mdct->length);
    zero_float(tmp + mdct->length * 2 - mdct->leading_zeros, mdct->leading_zeros);
    move_float(mdct->mem, tmp + mdct->length, mdct->mem_length);

    mult_vec(tmp, mdct->window, mdct->length * 2);

    LC3_INT hlen = mdct->length / 2;
    for (i = 0; i < hlen; i++) {
        output[i]        = -tmp[hlen * 3 - i - 1] - tmp[hlen * 3 + i];
        output[hlen + i] = tmp[i] - tmp[hlen * 2 - i - 1];
    }

    move_float(tmp, output, mdct->length);

    dct4_apply(&mdct->dct, tmp, output);
}

void processMdct_fl(LC3_FLOAT* in, LC3_FLOAT* out, Mdct* mdctStruct) { mdct_apply(in, out, mdctStruct); }
