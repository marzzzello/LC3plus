/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processResidualDecoding_fl(LC3_INT* bitsRead, LC3_FLOAT x[], LC3_INT L_spec, uint8_t prm[], LC3_INT resQBits
#ifdef ENABLE_HR_MODE
		, LC3_INT hrmode
#endif
)
{
    LC3_INT k = 0, n = 0;
    LC3_FLOAT offset1 = 0, offset2 = 0;
#ifdef ENABLE_HR_MODE
    LC3_FLOAT offset = 0;
#endif

#ifdef ENABLE_HR_MODE
    LC3_INT iter = 0, iter_max = 1;

    if (hrmode)
    {
    	iter_max = EXT_RES_ITER_MAX;
    	offset = offset1 = offset2 = 0.25;
    }
    else
    {
    	offset1 = 0.1875;
    	offset2 = 0.3125;
    }
#else
    offset1 = 0.1875;
    offset2 = 0.3125;
#endif

#ifdef ENABLE_HR_MODE
	if (hrmode)
	{
		while (n < resQBits && iter < iter_max)
		{
		    k = 0;
			while (k < L_spec && n < resQBits)
			{
				if (k < L_spec)
				{
					if (x[k] != 0)
					{
						if ((prm[n >> 3] & 1 << (n & 7)) == 0)
						{
							x[k] -= offset;
						}
						else
						{

							x[k] += offset;
						}
						if (++n >= resQBits)
						{
							break;
						}
					}
				}

				k ++;
			}
			offset /= 2;
			iter ++;
		}
	}
	else
	{
#endif
		while (k < L_spec && n < resQBits) {
			if (x[k] != 0) {
				if ((prm[n >> 3] & 1 << (n & 7)) == 0)
				{
					if (x[k] > 0) {
						x[k] -= offset1;
					} else {
						x[k] -= offset2;
					}
				} else {
					if (x[k] > 0) {
						x[k] += offset2;
					} else {
						x[k] += offset1;
					}
				}
				n++;
			}

			k++;
		}
#ifdef ENABLE_HR_MODE
	}
#endif
    *bitsRead = n;
}
