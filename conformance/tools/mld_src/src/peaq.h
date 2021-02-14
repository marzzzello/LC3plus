/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef PEAQ_H
#define PEAQ_H

#include <stdio.h>

#define PEAQ_SAMPLERATE 48000 // expected input sapmling rate

typedef struct peaq Peaq;

// allocate peaq struct, for normal use cases, playback level should be 92
Peaq* peaq_init(double playback_level);

// free peaq struct
void peaq_free(Peaq* pq);

// update peaq with pcm data. The is expected to be in ragne [-1,1]
void peaq_update(Peaq* pq, const float* pcm, int len);

// call finish once after last update to flush internal buffers
void peaq_finish(Peaq* pq);

// print maximum loudness difference between to files
// segment size is the number of frames that should be combined
void  peaq_print_mld(Peaq* pq1, Peaq* pq2, int segment_size, FILE* segment_outfile);
float peaq_get_mld(Peaq* pq1, Peaq* pq2);

#endif
