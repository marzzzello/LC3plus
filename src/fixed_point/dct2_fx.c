/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


void dct16_fx(const Word16 *in, Word16 *out)
{
    Dyn_Mem_Deluxe_In(
        Word16 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
        Word16 b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;
    );

    a0  = add(in[15], in[0]);
    a1  = add(in[14], in[1]);
    a2  = add(in[13], in[2]);
    a3  = add(in[12], in[3]);
    a4  = add(in[11], in[4]);
    a5  = add(in[10], in[5]);
    a6  = add(in[9], in[6]);
    a7  = add(in[8], in[7]);
    a10 = sub(in[5], in[10]);
    a11 = sub(in[4], in[11]);
    a12 = sub(in[3], in[12]);
    a13 = sub(in[2], in[13]);

    b0  = add(a7, a0);
    b1  = add(a6, a1);
    b2  = add(a5, a2);
    b3  = add(a4, a3);
    b4  = sub(a3, a4);
    b5  = sub(a2, a5);
    b6  = sub(a1, a6);
    b7  = sub(a0, a7);
    b8  = sub(in[7], in[8]);
    b9  = sub(in[6], in[9]);
    b10 = add(mult_r(a10, -23170), mult_r(a13, 23170)); /* -Cπ/4 Cπ/4 */
    b11 = add(mult_r(a11, -23170), mult_r(a12, 23170)); /* -Cπ/4 Cπ/4 */
    b12 = add(mult_r(a12, 23170), mult_r(a11, 23170));  /*  Cπ/4 Cπ/4 */
    b13 = add(mult_r(a13, 23170), mult_r(a10, 23170));  /*  Cπ/4 Cπ/4 */
    b14 = sub(in[1], in[14]);
    b15 = sub(in[0], in[15]);

    a0  = add(b3, b0);
    a1  = add(b2, b1);
    a2  = sub(b1, b2);
    a3  = sub(b0, b3);
    a4  = b4;                                         move16();
    a5  = add(mult_r(b5, -23170), mult_r(b6, 23170)); /* -Cπ/4 Cπ/4 */
    a6  = add(mult_r(b6, 23170), mult_r(b5, 23170));  /*  Cπ/4 Cπ/4 */
    a7  = b7;                                         move16();
    a8  = add(b11, b8);
    a9  = add(b10, b9);
    a10 = sub(b9, b10);
    a11 = sub(b8, b11);
    a12 = sub(b15, b12);
    a13 = sub(b14, b13);
    a14 = add(b13, b14);
    a15 = add(b12, b15);

    out[0]  = add(mult_r(a0, 8192), mult_r(a1, 8192)); move16();   /*  Cπ/4/√8   Cπ/4/√8  */
    out[8]  = add(mult_r(a1, -8192), mult_r(a0, 8192)); move16();  /* -Cπ/4/√8   Cπ/4/√8  */
    out[4]  = add(mult_r(a2, 4433), mult_r(a3, 10703)); move16();  /*  Sπ/8/√8   Cπ/8/√8  */
    out[12] = add(mult_r(a3, 4433), mult_r(a2, -10703)); move16(); /*  C3π/8/√8 -S3π/8/√8 */
    b4      = add(a5, a4);
    b5      = sub(a4, a5);
    b6      = sub(a7, a6);
    b7      = add(a6, a7);
    b8      = a8;                                            move16();
    b9      = add(mult_r(a9, -30274), mult_r(a14, 12540));   /* -Cπ/8  Sπ/8 */
    b10     = add(mult_r(a10, -12540), mult_r(a13, -30274)); /* -Sπ/8 -Cπ/8 */
    b11     = a11;                                           move16();
    b12     = a12;                                           move16();
    b13     = add(mult_r(a13, 12540), mult_r(a10, -30274));  /* C3π/8 -S3π/8 */
    b14     = add(mult_r(a14, 30274), mult_r(a9, 12540));    /* S3π/8  C3π/8 */
    b15     = a15;                                           move16();

    out[2]  = add(mult_r(b4, 2260), mult_r(b7, 11363)); move16();  /* Sπ/16/√8   Cπ/16/√8  */
    out[10] = add(mult_r(b5, 9633), mult_r(b6, 6436)); move16();   /* S5π/16/√8  C5π/16/√8 */
    out[6]  = add(mult_r(b6, 9633), mult_r(b5, -6436)); move16();  /* C3π/16/√8 -S3π/16/√8 */
    out[14] = add(mult_r(b7, 2260), mult_r(b4, -11363)); move16(); /* C7π/16/√8 -S7π/16/√8 */
    a8      = add(b9, b8);
    a9      = sub(b8, b9);
    a10     = sub(b11, b10);
    a11     = add(b10, b11);
    a12     = add(b13, b12);
    a13     = sub(b12, b13);
    a14     = sub(b15, b14);
    a15     = add(b14, b15);

    out[1]  = add(mult_r(a8, 1136), mult_r(a15, 11529)); move16();   /* Sπ/32/√8    Cπ/32/√8   */
    out[9]  = add(mult_r(a9, 8956), mult_r(a14, 7350)); move16();    /* S9π/32/√8   C9π/32/√8  */
    out[5]  = add(mult_r(a10, 5461), mult_r(a13, 10217)); move16();  /* S5π/32/√8   C5π/32/√8  */
    out[13] = add(mult_r(a11, 11086), mult_r(a12, 3363)); move16();  /* S13π/32/√8  C13π/32/√8 */
    out[3]  = add(mult_r(a12, 11086), mult_r(a11, -3363)); move16(); /* C3π/32/√8  -S3π/32/√8  */
    out[11] = add(mult_r(a13, 5461), mult_r(a10, -10217)); move16(); /* C11π/32/√8 -S11π/32/√8 */
    out[7]  = add(mult_r(a14, 8956), mult_r(a9, -7350)); move16();   /* C7π/32/√8  -S7π/32/√8  */
    out[15] = add(mult_r(a15, 1136), mult_r(a8, -11529)); move16();  /* C15π/32/√8 -S15/32/√8  */

    Dyn_Mem_Deluxe_Out();
}

void idct16_fx(const Word16 *in, Word16 *out)
{
    Dyn_Mem_Deluxe_In(
        Word16 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
        Word16 b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;
    );

    a8  = add(mult_r(in[1], 1136), mult_r(in[15], -11529)); /* Sπ/32/√8   -S15π/32/√8 */
    a9  = add(mult_r(in[9], 8956), mult_r(in[7], -7350));   /* S9π/32/√8  -S7π/32/√8  */
    a10 = add(mult_r(in[5], 5461), mult_r(in[11], -10217)); /* S5π/32/√8  -S11π/32/√8 */
    a11 = add(mult_r(in[13], 11086), mult_r(in[3], -3363)); /* S13π/32/√8 -S3π/32/√8  */
    a12 = add(mult_r(in[3], 11086), mult_r(in[13], 3363));  /* C3π/32/√8   C13π/32/√8 */
    a13 = add(mult_r(in[11], 5461), mult_r(in[5], 10217));  /* C11π/32/√8  C5π/32/√8  */
    a14 = add(mult_r(in[7], 8956), mult_r(in[9], 7350));    /* C7π/32/√8   C9π/32/√8  */
    a15 = add(mult_r(in[15], 1136), mult_r(in[1], 11529));  /* C15π/32/√8  Cπ/32/√8   */

    b4  = add(mult_r(in[2], 2260), mult_r(in[14], -11363)); /* Sπ/16/√8  -S7π/16/√8 */
    b5  = add(mult_r(in[10], 9633), mult_r(in[6], -6436));  /* S5π/16/√8 -S3π/16/√8 */
    b6  = add(mult_r(in[6], 9633), mult_r(in[10], 6436));   /* C3π/16/√8  C5π/16/√8 */
    b7  = add(mult_r(in[14], 2260), mult_r(in[2], 11363));  /* C7π/16/√8  Cπ/16/√8  */
    b8  = add(a9, a8);
    b9  = sub(a8, a9);
    b10 = sub(a11, a10);
    b11 = add(a10, a11);
    b12 = add(a13, a12);
    b13 = sub(a12, a13);
    b14 = sub(a15, a14);
    b15 = add(a14, a15);

    a0  = add(mult_r(in[0], 8192), mult_r(in[8], 8192));    /*  Cπ/4/√8  Cπ/4/√8  */
    a1  = add(mult_r(in[8], -8192), mult_r(in[0], 8192));   /* -Cπ/4/√8  Cπ/4/√8  */
    a2  = add(mult_r(in[4], 4433), mult_r(in[12], -10703)); /*  Sπ/8/√8 -S3π/8/√8 */
    a3  = add(mult_r(in[12], 4433), mult_r(in[4], 10703));  /*  C3π/8/√8 Cπ/8/√8  */
    a4  = add(b5, b4);
    a5  = sub(b4, b5);
    a6  = sub(b7, b6);
    a7  = add(b6, b7);
    a8  = b8;                                            move16();
    a9  = add(mult_r(b9, -30274), mult_r(b14, 12540));   /* -Cπ/8  C3π/8 */
    a10 = add(mult_r(b10, -12540), mult_r(b13, -30274)); /* -Sπ/8 -S3π/8 */
    a11 = b11;                                           move16();
    a12 = b12;                                           move16();
    a13 = add(mult_r(b13, 12540), mult_r(b10, -30274));  /* C3π/8 -Cπ/8 */
    a14 = add(mult_r(b14, 30274), mult_r(b9, 12540));    /* S3π/8  Sπ/8 */
    a15 = b15;                                           move16();

    b0  = add(a3, a0);
    b1  = add(a2, a1);
    b2  = sub(a1, a2);
    b3  = sub(a0, a3);
    b4  = a4;                                         move16();
    b5  = add(mult_r(a5, -23170), mult_r(a6, 23170)); /* -Cπ/4 Cπ/4 */
    b6  = add(mult_r(a6, 23170), mult_r(a5, 23170));  /*  Cπ/4 Cπ/4 */
    b7  = a7;                                         move16();
    b8  = add(a11, a8);
    b9  = add(a10, a9);
    b10 = sub(a9, a10);
    b11 = sub(a8, a11);
    b12 = sub(a15, a12);
    b13 = sub(a14, a13);
    b14 = add(a13, a14);
    b15 = add(a12, a15);

    a0  = add(b7, b0);
    a1  = add(b6, b1);
    a2  = add(b5, b2);
    a3  = add(b4, b3);
    a4  = sub(b3, b4);
    a5  = sub(b2, b5);
    a6  = sub(b1, b6);
    a7  = sub(b0, b7);
    a10 = add(mult_r(b10, -23170), mult_r(b13, 23170)); /* -Cπ/4 Cπ/4 */
    a11 = add(mult_r(b11, -23170), mult_r(b12, 23170)); /* -Cπ/4 Cπ/4 */
    a12 = add(mult_r(b12, 23170), mult_r(b11, 23170));  /*  Cπ/4 Cπ/4 */
    a13 = add(mult_r(b13, 23170), mult_r(b10, 23170));  /*  Cπ/4 Cπ/4 */

    out[0]  = add(b15, a0); move16();
    out[1]  = add(b14, a1); move16();
    out[2]  = add(a13, a2); move16();
    out[3]  = add(a12, a3); move16();
    out[4]  = add(a11, a4); move16();
    out[5]  = add(a10, a5); move16();
    out[6]  = add(b9, a6);  move16();
    out[7]  = add(b8, a7);  move16();
    out[8]  = sub(a7, b8);  move16();
    out[9]  = sub(a6, b9);  move16();
    out[10] = sub(a5, a10); move16();
    out[11] = sub(a4, a11); move16();
    out[12] = sub(a3, a12); move16();
    out[13] = sub(a2, a13); move16();
    out[14] = sub(a1, b14); move16();
    out[15] = sub(a0, b15); move16();

    Dyn_Mem_Deluxe_Out();
}

