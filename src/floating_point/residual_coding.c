/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processResidualCoding_fl(LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gain, LC3_INT L_spec, LC3_INT targetBits, LC3_INT nBits, uint8_t* resBits, LC3_INT* numResBits
#ifdef ENABLE_HR_MODE
		, LC3_INT hrmode
#endif
)
{
    LC3_INT n = 0, m = 0, k = 0;
#ifdef ENABLE_HR_MODE
    LC3_INT iter=0;
    LC3_FLOAT offset;
    LC3_INT iter_max = 1;
#endif

#ifdef ENABLE_HR_MODE

    memset(resBits, 0, MAX_RESBITS_LEN);

    m = targetBits - nBits + 4;
    if (hrmode)
    {
        m += 10;
    }

    assert(m <= MAX_RESBITS);

    offset = .25;
    if (hrmode)
    {
    	iter_max = EXT_RES_ITER_MAX;
    }
    while (iter < iter_max && n < m)
    {
    	k = 0;
		while (k < L_spec && n < m)
		{
			if (xq[k] != 0)
			{
				if (x[k] >= (LC3_FLOAT)xq[k] * gain)
				{
					resBits[n >> 3] |= 1 << (n & 7);
					x[k] -= gain * offset;
				}
				else
				{
					resBits[n >> 3] &= ~(1 << (n & 7));
					x[k] += gain * offset;
				}

				n++;
			}

			k++;
		}
	iter ++;
	offset *= .5;
    }
#else
    m = targetBits - nBits + 4;

    while (k < L_spec && n < m) {
        if (xq[k] != 0) {
            if (x[k] >= (LC3_FLOAT)xq[k] * gain) {
            	resBits[n >> 3] |= 1 << (n & 7);
            } else {
            	resBits[n >> 3] &= ~(1 << (n & 7));
            }

            n++;
        }

        k++;
    }
#endif

    *numResBits = n;
}
