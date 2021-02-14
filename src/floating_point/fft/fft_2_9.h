/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

/* guard against unindended includes */
#ifndef INCLUDED_FROM_IISFFT_C
#error "this file must not be included"
#endif

static __inline void butterfly(const LC3_FLOAT a, const LC3_FLOAT b, LC3_FLOAT* aPlusb, LC3_FLOAT* aMinusb)
{
    *aPlusb = a + b;
    *aMinusb = a - b;
}

static void fft2(LC3_FLOAT* vec)
{
    /*
    1.0000  1.0000
    1.0000 -1.0000
    */
    LC3_FLOAT re1 = vec[0];
    LC3_FLOAT im1 = vec[1];
    LC3_FLOAT re2 = vec[2];
    LC3_FLOAT im2 = vec[3];

    vec[0] = re1 + re2;
    vec[1] = im1 + im2;
    vec[2] = re1 - re2;
    vec[3] = im1 - im2;
}

static void fft3(LC3_FLOAT* vec)
{
    const LC3_FLOAT C31 = 0.5;               /* cos(PI/3); sin(2*PI/3) */
    const LC3_FLOAT C32 = 0.866025403784439; /* cos(PI/3); sin(2*PI/3) */

    LC3_FLOAT re1 = vec[0];
    LC3_FLOAT im1 = vec[1];
    LC3_FLOAT re2 = vec[2];
    LC3_FLOAT im2 = vec[3];
    LC3_FLOAT re3 = vec[4];
    LC3_FLOAT im3 = vec[5];
    /*
    1.0000             1.0000             1.0000
                       C31                C32
    1.0000            -0.5000 - 0.8660i  -0.5000 + 0.8660i
    1.0000            -0.5000 + 0.8660i  -0.5000 - 0.8660i
    */
    LC3_FLOAT tmp1 = re2 + re3;
    LC3_FLOAT tmp3 = im2 + im3;
    LC3_FLOAT tmp2 = re2 - re3;
    LC3_FLOAT tmp4 = im2 - im3;

    vec[0] = re1 + tmp1;
    vec[1] = im1 + tmp3;
    vec[2] = re1 - C31 * tmp1 + C32 * tmp4;
    vec[4] = re1 - C31 * tmp1 - C32 * tmp4;
    vec[3] = im1 - C32 * tmp2 - C31 * tmp3;
    vec[5] = im1 + C32 * tmp2 - C31 * tmp3;
}

static void fft4(LC3_FLOAT* vec)
{
    LC3_FLOAT temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;

    /* Pre-additions */
    temp0 = vec[0] + vec[4];
    temp2 = vec[0] - vec[4];
    temp1 = vec[1] + vec[5];
    temp3 = vec[1] - vec[5];
    temp4 = vec[2] + vec[6];
    temp7 = vec[2] - vec[6];
    temp5 = vec[7] + vec[3];
    temp6 = vec[7] - vec[3];

    /* Post-additions */
    vec[0] = temp0 + temp4;
    vec[1] = temp1 + temp5;
    vec[2] = temp2 - temp6;
    vec[3] = temp3 - temp7;
    vec[4] = temp0 - temp4;
    vec[5] = temp1 - temp5;
    vec[6] = temp2 + temp6;
    vec[7] = temp3 + temp7;
}

static void fft5(LC3_FLOAT* vec)
{
    const LC3_FLOAT C51 = 0.309016994374947; /* cos(2*PI/5);   */
    const LC3_FLOAT C52 = 0.951056516295154; /* sin(2*PI/5);   */
    const LC3_FLOAT C53 = 0.809016994374947; /* cos(  PI/5);   */
    const LC3_FLOAT C54 = 0.587785252292473; /* sin(  PI/5);   */

    LC3_FLOAT re1, im1, re2, im2, re3, im3, re4, im4, re5, im5, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

    re1 = vec[0];
    im1 = vec[1];
    re2 = vec[2];
    im2 = vec[3];
    re3 = vec[4];
    im3 = vec[5];
    re4 = vec[6];
    im4 = vec[7];
    re5 = vec[8];
    im5 = vec[9];
    /*
    1.0000             1.0000             1.0000             1.0000             1.0000
                       C51                C52                C53                C54
    1.0000             0.3090 - 0.9511i  -0.8090 - 0.5878i  -0.8090 + 0.5878i   0.3090 + 0.9511i
    1.0000            -0.8090 - 0.5878i   0.3090 + 0.9511i   0.3090 - 0.9511i  -0.8090 + 0.5878i
    1.0000            -0.8090 + 0.5878i   0.3090 - 0.9511i   0.3090 + 0.9511i  -0.8090 - 0.5878i
    1.0000             0.3090 + 0.9511i  -0.8090 + 0.5878i  -0.8090 - 0.5878i   0.3090 - 0.9511i
    */
    tmp1 = re2 + re5;
    tmp2 = re2 - re5;
    tmp3 = im2 + im5;
    tmp4 = im2 - im5;
    tmp5 = re3 + re4;
    tmp6 = re3 - re4;
    tmp7 = im3 + im4;
    tmp8 = im3 - im4;

    vec[0] = re1 + tmp1 + tmp5;
    vec[1] = im1 + tmp3 + tmp7;
    vec[2] = re1 + C51 * tmp1 - C53 * tmp5 + C52 * tmp4 + C54 * tmp8;
    vec[8] = re1 + C51 * tmp1 - C53 * tmp5 - C52 * tmp4 - C54 * tmp8;
    vec[3] = im1 - C52 * tmp2 - C54 * tmp6 + C51 * tmp3 - C53 * tmp7;
    vec[9] = im1 + C52 * tmp2 + C54 * tmp6 + C51 * tmp3 - C53 * tmp7;
    vec[4] = re1 - C53 * tmp1 + C51 * tmp5 + C54 * tmp4 - C52 * tmp8;
    vec[6] = re1 - C53 * tmp1 + C51 * tmp5 - C54 * tmp4 + C52 * tmp8;
    vec[5] = im1 - C54 * tmp2 + C52 * tmp6 - C53 * tmp3 + C51 * tmp7;
    vec[7] = im1 + C54 * tmp2 - C52 * tmp6 - C53 * tmp3 + C51 * tmp7;
}

static void ifft5(LC3_FLOAT* vec)
{
    const LC3_FLOAT C51 = 0.309016994374947; /* cos(2*PI/5);   */
    const LC3_FLOAT C52 = 0.951056516295154; /* sin(2*PI/5);   */
    const LC3_FLOAT C53 = 0.809016994374947; /* cos(  PI/5);   */
    const LC3_FLOAT C54 = 0.587785252292473; /* sin(  PI/5);   */

    LC3_FLOAT re1, im1, re2, im2, re3, im3, re4, im4, re5, im5, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

    re1 = vec[0];
    im1 = vec[1];
    re2 = vec[2];
    im2 = vec[3];
    re3 = vec[4];
    im3 = vec[5];
    re4 = vec[6];
    im4 = vec[7];
    re5 = vec[8];
    im5 = vec[9];
    /*
    1.0000             1.0000             1.0000             1.0000             1.0000
                       C51      C52       C53      C54
    1.0000             0.3090 + 0.9511i  -0.8090 + 0.5878i  -0.8090 - 0.5878i   0.3090 - 0.9511i
    1.0000            -0.8090 + 0.5878i   0.3090 - 0.9511i   0.3090 + 0.9511i  -0.8090 - 0.5878i
    1.0000            -0.8090 - 0.5878i   0.3090 + 0.9511i   0.3090 - 0.9511i  -0.8090 + 0.5878i
    1.0000             0.3090 - 0.9511i  -0.8090 - 0.5878i  -0.8090 + 0.5878i   0.3090 + 0.9511i
    */
    tmp1 = re2 + re5;
    tmp2 = re2 - re5;
    tmp3 = im2 + im5;
    tmp4 = im2 - im5;
    tmp5 = re3 + re4;
    tmp6 = re3 - re4;
    tmp7 = im3 + im4;
    tmp8 = im3 - im4;

    vec[0] = re1 + tmp1 + tmp5;
    vec[1] = im1 + tmp3 + tmp7;
    vec[2] = re1 + C51 * tmp1 - C53 * tmp5 - C52 * tmp4 - C54 * tmp8;
    vec[3] = im1 + C52 * tmp2 + C54 * tmp6 + C51 * tmp3 - C53 * tmp7;
    vec[4] = re1 - C53 * tmp1 + C51 * tmp5 - C54 * tmp4 + C52 * tmp8;
    vec[5] = im1 + C54 * tmp2 - C52 * tmp6 - C53 * tmp3 + C51 * tmp7;
    vec[6] = re1 - C53 * tmp1 + C51 * tmp5 + C54 * tmp4 - C52 * tmp8;
    vec[7] = im1 - C54 * tmp2 + C52 * tmp6 - C53 * tmp3 + C51 * tmp7;
    vec[8] = re1 + C51 * tmp1 - C53 * tmp5 + C52 * tmp4 + C54 * tmp8;
    vec[9] = im1 - C52 * tmp2 - C54 * tmp6 + C51 * tmp3 - C53 * tmp7;
}

static void fft7(LC3_FLOAT* vec)
{
    const LC3_FLOAT C71R = 0.623489801858734;
    const LC3_FLOAT C72R = -0.222520933956314;
    const LC3_FLOAT C73R = -0.900968867902419;
    const LC3_FLOAT C71I = -0.78183148246803;
    const LC3_FLOAT C72I = -0.974927912181824;
    const LC3_FLOAT C73I = -0.433883739117558;

    LC3_FLOAT re0, im0, re1_6p, re1_6m, im1_6p, im1_6m, re2_5p, re2_5m, im2_5p, im2_5m, re3_4p, re3_4m, im3_4p, im3_4m;

    re0 = vec[0];
    im0 = vec[1];

    butterfly(vec[1 * 2], vec[6 * 2], &re1_6p, &re1_6m);
    butterfly(vec[1 * 2 + 1], vec[6 * 2 + 1], &im1_6p, &im1_6m);
    butterfly(vec[2 * 2], vec[5 * 2], &re2_5p, &re2_5m);
    butterfly(vec[2 * 2 + 1], vec[5 * 2 + 1], &im2_5p, &im2_5m);
    butterfly(vec[3 * 2], vec[4 * 2], &re3_4p, &re3_4m);
    butterfly(vec[3 * 2 + 1], vec[4 * 2 + 1], &im3_4p, &im3_4m);

    vec[0] = re0 + re1_6p + re2_5p + re3_4p;
    vec[1] = im0 + im1_6p + im2_5p + im3_4p;
    /*
        1 + 0i   1 + 0i          1 + 0i          1 + 0i          1 + 0i          1 + 0i	        1 + 0i
        1 + 0i	 C71R - C71Ii	-C72R - C72Ii	-C73R - C73Ii	-C73R + C73Ii	-C72R + C72Ii	C71R + C71Ii
        1 + 0i	 C72R - C72Ii	-C73R + C73Ii	 C71R + C71Ii    C71R - C71Ii	-C73R - C73Ii  -C72R + C72Ii
        1 + 0i	 C73R - C73Ii	 C71R + C71Ii	-C72R - C72Ii	-C72R + C72Ii	 C71R - C71Ii  -C73R + C73Ii
        1 + 0i	 C73R + C73Ii	 C71R - C71Ii	-C72R + C72Ii	-C72R - C72Ii	 C71R + C71Ii  -C73R - C73Ii
        1 + 0i	 C72R + C72Ii	-C73R - C73Ii	 C71R - C71Ii	 C71R + C71Ii	-C73R + C73Ii  -C72R - C72Ii
        1 + 0i	 C71R + C71Ii	-C72R + C72Ii	-C73R + C73Ii	-C73R - C73Ii	-C72R - C72Ii   C71R - C71Ii
    */
    vec[2] = re0 + C71R * re1_6p + C72R * re2_5p + C73R * re3_4p - C71I * im1_6m - C72I * im2_5m - C73I * im3_4m;
    vec[12] = re0 + C71R * re1_6p + C72R * re2_5p + C73R * re3_4p + C71I * im1_6m + C72I * im2_5m + C73I * im3_4m;
    vec[4] = re0 + C72R * re1_6p + C73R * re2_5p + C71R * re3_4p - C72I * im1_6m + C73I * im2_5m + C71I * im3_4m;
    vec[10] = re0 + C72R * re1_6p + C73R * re2_5p + C71R * re3_4p + C72I * im1_6m - C73I * im2_5m - C71I * im3_4m;
    vec[6] = re0 + C73R * re1_6p + C71R * re2_5p + C72R * re3_4p - C73I * im1_6m + C71I * im2_5m - C72I * im3_4m;
    vec[8] = re0 + C73R * re1_6p + C71R * re2_5p + C72R * re3_4p + C73I * im1_6m - C71I * im2_5m + C72I * im3_4m;

    vec[3] = im0 + C71I * re1_6m + C72I * re2_5m + C73I * re3_4m + C71R * im1_6p + C72R * im2_5p + C73R * im3_4p;
    vec[13] = im0 - C71I * re1_6m - C72I * re2_5m - C73I * re3_4m + C71R * im1_6p + C72R * im2_5p + C73R * im3_4p;
    vec[5] = im0 + C72I * re1_6m - C73I * re2_5m - C71I * re3_4m + C72R * im1_6p + C73R * im2_5p + C71R * im3_4p;
    vec[11] = im0 - C72I * re1_6m + C73I * re2_5m + C71I * re3_4m + C72R * im1_6p + C73R * im2_5p + C71R * im3_4p;
    vec[7] = im0 + C73I * re1_6m - C71I * re2_5m + C72I * re3_4m + C73R * im1_6p + C71R * im2_5p + C72R * im3_4p;
    vec[9] = im0 - C73I * re1_6m + C71I * re2_5m - C72I * re3_4m + C73R * im1_6p + C71R * im2_5p + C72R * im3_4p;
}

static void fft8(LC3_FLOAT* vec)
{
    const LC3_FLOAT INV_SQRT2 = 7.071067811865475e-1;
    LC3_FLOAT temp1[16], temp2[16];

    /* Pre-additions */
    temp1[0] = vec[0] + vec[8];
    temp1[2] = vec[0] - vec[8];
    temp1[1] = vec[1] + vec[9];
    temp1[3] = vec[1] - vec[9];
    temp1[4] = vec[2] + vec[10];
    temp1[6] = vec[2] - vec[10];
    temp1[5] = vec[3] + vec[11];
    temp1[7] = vec[3] - vec[11];
    temp1[8] = vec[4] + vec[12];
    temp1[10] = vec[4] - vec[12];
    temp1[9] = vec[5] + vec[13];
    temp1[11] = vec[5] - vec[13];
    temp1[12] = vec[6] + vec[14];
    temp1[14] = vec[6] - vec[14];
    temp1[13] = vec[7] + vec[15];
    temp1[15] = vec[7] - vec[15];

    /* Pre-additions and core multiplications */
    temp2[0] = temp1[0] + temp1[8];
    temp2[4] = temp1[0] - temp1[8];
    temp2[1] = temp1[1] + temp1[9];
    temp2[5] = temp1[1] - temp1[9];
    temp2[8] = temp1[2] - temp1[11];
    temp2[10] = temp1[2] + temp1[11];
    temp2[9] = temp1[3] + temp1[10];
    temp2[11] = temp1[3] - temp1[10];
    temp2[2] = temp1[4] + temp1[12];
    temp2[7] = temp1[4] - temp1[12];
    temp2[3] = temp1[5] + temp1[13];
    temp2[6] = temp1[13] - temp1[5];

    temp1[1] = temp1[6] + temp1[14];
    temp1[2] = temp1[6] - temp1[14];
    temp1[0] = temp1[7] + temp1[15];
    temp1[3] = temp1[7] - temp1[15];
    temp2[12] = (temp1[0] + temp1[2]) * INV_SQRT2;
    temp2[14] = (temp1[0] - temp1[2]) * INV_SQRT2;
    temp2[13] = (temp1[3] - temp1[1]) * INV_SQRT2;
    temp2[15] = (temp1[1] + temp1[3]) * -INV_SQRT2;

    /* Post-additions */
    vec[0] = temp2[0] + temp2[2];
    vec[8] = temp2[0] - temp2[2];
    vec[1] = temp2[1] + temp2[3];
    vec[9] = temp2[1] - temp2[3];
    vec[4] = temp2[4] - temp2[6];
    vec[12] = temp2[4] + temp2[6];
    vec[5] = temp2[5] - temp2[7];
    vec[13] = temp2[5] + temp2[7];
    vec[6] = temp2[8] + temp2[14];
    vec[14] = temp2[8] - temp2[14];
    vec[7] = temp2[9] + temp2[15];
    vec[15] = temp2[9] - temp2[15];
    vec[2] = temp2[10] + temp2[12];
    vec[10] = temp2[10] - temp2[12];
    vec[3] = temp2[11] + temp2[13];
    vec[11] = temp2[11] - temp2[13];
}

static void ifft8(LC3_FLOAT* vec)
{
    const LC3_FLOAT INV_SQRT2 = 7.071067811865475e-1;

    LC3_FLOAT temp1[16], temp2[16];

    /* Pre-additions */
    temp1[0] = vec[0] + vec[8];
    temp1[2] = vec[0] - vec[8];
    temp1[1] = vec[1] + vec[9];
    temp1[3] = vec[1] - vec[9];
    temp1[4] = vec[2] + vec[10];
    temp1[6] = vec[2] - vec[10];
    temp1[5] = vec[3] + vec[11];
    temp1[7] = vec[3] - vec[11];
    temp1[8] = vec[4] + vec[12];
    temp1[10] = vec[4] - vec[12];
    temp1[9] = vec[5] + vec[13];
    temp1[11] = vec[5] - vec[13];
    temp1[12] = vec[6] + vec[14];
    temp1[14] = vec[6] - vec[14];
    temp1[13] = vec[7] + vec[15];
    temp1[15] = vec[7] - vec[15];

    /* Pre-additions and core multiplications */
    temp2[0] = temp1[0] + temp1[8];
    temp2[4] = temp1[0] - temp1[8];
    temp2[1] = temp1[1] + temp1[9];
    temp2[5] = temp1[1] - temp1[9];
    temp2[8] = temp1[2] - temp1[11];
    temp2[10] = temp1[2] + temp1[11];
    temp2[9] = temp1[3] + temp1[10];
    temp2[11] = temp1[3] - temp1[10];
    temp2[2] = temp1[4] + temp1[12];
    temp2[7] = temp1[4] - temp1[12];
    temp2[3] = temp1[5] + temp1[13];
    temp2[6] = temp1[13] - temp1[5];

    temp1[1] = temp1[6] + temp1[14];
    temp1[2] = temp1[6] - temp1[14];
    temp1[0] = temp1[7] + temp1[15];
    temp1[3] = temp1[7] - temp1[15];
    temp2[12] = (temp1[0] + temp1[2]) * INV_SQRT2;
    temp2[14] = (temp1[0] - temp1[2]) * INV_SQRT2;
    temp2[13] = (temp1[3] - temp1[1]) * INV_SQRT2;
    temp2[15] = (temp1[1] + temp1[3]) * -INV_SQRT2;

    /* Post-additions */
    vec[0] = temp2[0] + temp2[2];
    vec[8] = temp2[0] - temp2[2];
    vec[1] = temp2[1] + temp2[3];
    vec[9] = temp2[1] - temp2[3];
    vec[4] = temp2[4] + temp2[6];
    vec[12] = temp2[4] - temp2[6];
    vec[5] = temp2[5] + temp2[7];
    vec[13] = temp2[5] - temp2[7];
    vec[6] = temp2[10] - temp2[12];
    vec[14] = temp2[10] + temp2[12];
    vec[7] = temp2[11] - temp2[13];
    vec[15] = temp2[11] + temp2[13];
    vec[2] = temp2[8] - temp2[14];
    vec[10] = temp2[8] + temp2[14];
    vec[3] = temp2[9] - temp2[15];
    vec[11] = temp2[9] + temp2[15];
}

static void fft9(LC3_FLOAT* vec)
{
    const LC3_FLOAT C91 = 0.766044443118978; /* cos(2*PI/5);   */
    const LC3_FLOAT C92 = 0.642787609686539;
    const LC3_FLOAT C93 = 0.17364817766693;
    const LC3_FLOAT C94 = 0.984807753012208;
    const LC3_FLOAT C95 = 0.5;
    const LC3_FLOAT C96 = 0.866025403784439;
    const LC3_FLOAT C97 = 0.939692620785908;
    const LC3_FLOAT C98 = 0.342020143325669;

    LC3_FLOAT re1, im1, re2_9p, re2_9m, im2_9p, im2_9m, re3_8p, re3_8m, im3_8p, im3_8m, re4_7p, re4_7m, im4_7p, im4_7m,
        re5_6p, re5_6m, im5_6p, im5_6m;

    re1 = vec[0];
    im1 = vec[1];

    butterfly(vec[1 * 2], vec[8 * 2], &re2_9p, &re2_9m);
    butterfly(vec[1 * 2 + 1], vec[8 * 2 + 1], &im2_9p, &im2_9m);
    butterfly(vec[2 * 2], vec[7 * 2], &re3_8p, &re3_8m);
    butterfly(vec[2 * 2 + 1], vec[7 * 2 + 1], &im3_8p, &im3_8m);
    butterfly(vec[3 * 2], vec[6 * 2], &re4_7p, &re4_7m);
    butterfly(vec[3 * 2 + 1], vec[6 * 2 + 1], &im4_7p, &im4_7m);
    butterfly(vec[4 * 2], vec[5 * 2], &re5_6p, &re5_6m);
    butterfly(vec[4 * 2 + 1], vec[5 * 2 + 1], &im5_6p, &im5_6m);

    /*
    1.0000     1.0000       1.0000       1.0000       1.0000       1.0000       1.0000       1.0000       1.0000
    1.0000   C91 - C92i   C93 - C94i  -C95 - C96i  -C97 - C98i  -C97 + C98i  -C95 + C96i   C93 + C94i   C91 + C92i
    1.0000   C93 - C94i  -C97 - C98i  -C95 + C96i   C91 + C92i   C91 - C92i  -C95 - C96i  -C97 + C98i   C93 + C94i
    1.0000  -C95 - C96i  -C95 + C96i     1.0000    -C95 - C96i  -C95 + C96i     1.0000    -C95 - C96i  -C95 + C96i
    1.0000  -C97 - C98i   C91 + C92i  -C95 - C96i   C93 + C94i   C93 - C94i  -C95 + C96i   C91 - C92i  -C97 + C98i
    1.0000  -C97 + C98i   C91 - C92i  -C95 + C96i   C93 - C94i   C93 + C94i  -C95 - C96i   C91 + C92i  -C97 - C98i
    1.0000  -C95 + C96i  -C95 - C96i     1.0000    -C95 + C96i  -C95 - C96i     1.0000    -C95 + C96i  -C95 - C96i
    1.0000   C93 + C94i  -C97 + C98i  -C95 - C96i   C91 - C92i   C91 + C92i  -C95 + C96i  -C97 - C98i   C93 - C94i
    1.0000   C91 + C92i   C93 + C94i  -C95 + C96i  -C97 + C98i  -C97 - C98i  -C95 - C96i   C93 - C94i   C91 - C92i
    */
    vec[0] = re1 + re2_9p + re3_8p + re4_7p + re5_6p;
    vec[1] = im1 + im2_9p + im3_8p + im4_7p + im5_6p;
    vec[2] = re1 + C91 * re2_9p + C93 * re3_8p - C95 * re4_7p - C97 * re5_6p + C92 * im2_9m + C94 * im3_8m +
             C96 * im4_7m + C98 * im5_6m;
    vec[16] = re1 + C91 * re2_9p + C93 * re3_8p - C95 * re4_7p - C97 * re5_6p - C92 * im2_9m - C94 * im3_8m -
              C96 * im4_7m - C98 * im5_6m;
    vec[3] = im1 - C92 * re2_9m - C94 * re3_8m - C96 * re4_7m - C98 * re5_6m + C91 * im2_9p + C93 * im3_8p -
             C95 * im4_7p - C97 * im5_6p;
    vec[17] = im1 + C92 * re2_9m + C94 * re3_8m + C96 * re4_7m + C98 * re5_6m + C91 * im2_9p + C93 * im3_8p -
              C95 * im4_7p - C97 * im5_6p;
    vec[4] = re1 + C93 * re2_9p - C97 * re3_8p - C95 * re4_7p + C91 * re5_6p + C94 * im2_9m + C98 * im3_8m -
             C96 * im4_7m - C92 * im5_6m;
    vec[14] = re1 + C93 * re2_9p - C97 * re3_8p - C95 * re4_7p + C91 * re5_6p - C94 * im2_9m - C98 * im3_8m +
              C96 * im4_7m + C92 * im5_6m;
    vec[5] = im1 - C94 * re2_9m - C98 * re3_8m + C96 * re4_7m + C92 * re5_6m + C93 * im2_9p - C97 * im3_8p -
             C95 * im4_7p + C91 * im5_6p;
    vec[15] = im1 + C94 * re2_9m + C98 * re3_8m - C96 * re4_7m - C92 * re5_6m + C93 * im2_9p - C97 * im3_8p -
              C95 * im4_7p + C91 * im5_6p;
    vec[6] = re1 - C95 * (re2_9p + re3_8p + re5_6p) + re4_7p + C96 * (im2_9m - im3_8m + im5_6m);
    vec[12] = re1 - C95 * (re2_9p + re3_8p + re5_6p) + re4_7p - C96 * (im2_9m - im3_8m + im5_6m);
    vec[7] = im1 - C96 * (re2_9m - re3_8m + re5_6m) - C95 * (im2_9p + im3_8p + im5_6p) + im4_7p;
    vec[13] = im1 + C96 * (re2_9m - re3_8m + re5_6m) - C95 * (im2_9p + im3_8p + im5_6p) + im4_7p;
    vec[8] = re1 - C97 * re2_9p + C91 * re3_8p - C95 * re4_7p + C93 * re5_6p + C98 * im2_9m - C92 * im3_8m +
             C96 * im4_7m - C94 * im5_6m;
    vec[10] = re1 - C97 * re2_9p + C91 * re3_8p - C95 * re4_7p + C93 * re5_6p - C98 * im2_9m + C92 * im3_8m -
              C96 * im4_7m + C94 * im5_6m;
    vec[9] = im1 - C98 * re2_9m + C92 * re3_8m - C96 * re4_7m + C94 * re5_6m - C97 * im2_9p + C91 * im3_8p -
             C95 * im4_7p + C93 * im5_6p;
    vec[11] = im1 + C98 * re2_9m - C92 * re3_8m + C96 * re4_7m - C94 * re5_6m - C97 * im2_9p + C91 * im3_8p -
              C95 * im4_7p + C93 * im5_6p;
}

static void ifft9(LC3_FLOAT* vec)
{
    const LC3_FLOAT C91 = 0.766044443118978; /* cos(2*PI/5);   */
    const LC3_FLOAT C92 = 0.642787609686539;
    const LC3_FLOAT C93 = 0.17364817766693;
    const LC3_FLOAT C94 = 0.984807753012208;
    const LC3_FLOAT C95 = 0.5;
    const LC3_FLOAT C96 = 0.866025403784439;
    const LC3_FLOAT C97 = 0.939692620785908;
    const LC3_FLOAT C98 = 0.342020143325669;

    LC3_FLOAT re1, im1, re2_9p, re2_9m, im2_9p, im2_9m, re3_8p, re3_8m, im3_8p, im3_8m, re4_7p, re4_7m, im4_7p, im4_7m,
        re5_6p, re5_6m, im5_6p, im5_6m;

    re1 = vec[0];
    im1 = vec[1];
    butterfly(vec[1 * 2], vec[8 * 2], &re2_9p, &re2_9m);
    butterfly(vec[1 * 2 + 1], vec[8 * 2 + 1], &im2_9p, &im2_9m);
    butterfly(vec[2 * 2], vec[7 * 2], &re3_8p, &re3_8m);
    butterfly(vec[2 * 2 + 1], vec[7 * 2 + 1], &im3_8p, &im3_8m);
    butterfly(vec[3 * 2], vec[6 * 2], &re4_7p, &re4_7m);
    butterfly(vec[3 * 2 + 1], vec[6 * 2 + 1], &im4_7p, &im4_7m);
    butterfly(vec[4 * 2], vec[5 * 2], &re5_6p, &re5_6m);
    butterfly(vec[4 * 2 + 1], vec[5 * 2 + 1], &im5_6p, &im5_6m);

    /*
        1.0000     1.0000     1.0000       1.0000       1.0000       1.0000       1.0000       1.0000       1.0000
        1.0000   C91 - C92i   C93 - C94i  -C95 - C96i  -C97 - C98i  -C97 + C98i  -C95 + C96i   C93 + C94i   C91 + C92i
        1.0000   C93 - C94i  -C97 - C98i  -C95 + C96i   C91 + C92i   C91 - C92i  -C95 - C96i  -C97 + C98i   C93 + C94i
        1.0000  -C95 - C96i  -C95 + C96i     1.0000    -C95 - C96i  -C95 + C96i     1.0000    -C95 - C96i  -C95 + C96i
        1.0000  -C97 - C98i   C91 + C92i  -C95 - C96i   C93 + C94i   C93 - C94i  -C95 + C96i   C91 - C92i  -C97 + C98i
        1.0000  -C97 + C98i   C91 - C92i  -C95 + C96i   C93 - C94i   C93 + C94i  -C95 - C96i   C91 + C92i  -C97 - C98i
        1.0000  -C95 + C96i  -C95 - C96i     1.0000    -C95 + C96i  -C95 - C96i     1.0000    -C95 + C96i  -C95 - C96i
        1.0000   C93 + C94i  -C97 + C98i  -C95 - C96i   C91 - C92i   C91 + C92i  -C95 + C96i  -C97 - C98i   C93 - C94i
        1.0000   C91 + C92i   C93 + C94i  -C95 + C96i  -C97 + C98i  -C97 - C98i  -C95 - C96i   C93 - C94i   C91 - C92i
    */
    vec[0] = re1 + re2_9p + re3_8p + re4_7p + re5_6p;
    vec[1] = im1 + im2_9p + im3_8p + im4_7p + im5_6p;
    vec[2] = re1 + C91 * re2_9p + C93 * re3_8p - C95 * re4_7p - C97 * re5_6p - C92 * im2_9m - C94 * im3_8m -
             C96 * im4_7m - C98 * im5_6m;
    vec[3] = im1 + C92 * re2_9m + C94 * re3_8m + C96 * re4_7m + C98 * re5_6m + C91 * im2_9p + C93 * im3_8p -
             C95 * im4_7p - C97 * im5_6p;
    vec[4] = re1 + C93 * re2_9p - C97 * re3_8p - C95 * re4_7p + C91 * re5_6p - C94 * im2_9m - C98 * im3_8m +
             C96 * im4_7m + C92 * im5_6m;
    vec[5] = im1 + C94 * re2_9m + C98 * re3_8m - C96 * re4_7m - C92 * re5_6m + C93 * im2_9p - C97 * im3_8p -
             C95 * im4_7p + C91 * im5_6p;
    vec[6] = re1 - C95 * re2_9p - C95 * re3_8p + re4_7p - C95 * re5_6p - C96 * im2_9m + C96 * im3_8m - C96 * im5_6m;
    vec[7] = im1 + C96 * re2_9m - C96 * re3_8m + C96 * re5_6m - C95 * im2_9p - C95 * im3_8p + im4_7p - C95 * im5_6p;
    vec[8] = re1 - C97 * re2_9p + C91 * re3_8p - C95 * re4_7p + C93 * re5_6p - C98 * im2_9m + C92 * im3_8m -
             C96 * im4_7m + C94 * im5_6m;
    vec[9] = im1 + C98 * re2_9m - C92 * re3_8m + C96 * re4_7m - C94 * re5_6m - C97 * im2_9p + C91 * im3_8p -
             C95 * im4_7p + C93 * im5_6p;
    vec[10] = re1 - C97 * re2_9p + C91 * re3_8p - C95 * re4_7p + C93 * re5_6p + C98 * im2_9m - C92 * im3_8m +
              C96 * im4_7m - C94 * im5_6m;
    vec[11] = im1 - C98 * re2_9m + C92 * re3_8m - C96 * re4_7m + C94 * re5_6m - C97 * im2_9p + C91 * im3_8p -
              C95 * im4_7p + C93 * im5_6p;
    vec[12] = re1 - C95 * re2_9p - C95 * re3_8p + re4_7p - C95 * re5_6p + C96 * im2_9m - C96 * im3_8m + C96 * im5_6m;
    vec[13] = im1 - C96 * re2_9m + C96 * re3_8m - C96 * re5_6m - C95 * im2_9p - C95 * im3_8p + im4_7p - C95 * im5_6p;
    vec[14] = re1 + C93 * re2_9p - C97 * re3_8p - C95 * re4_7p + C91 * re5_6p + C94 * im2_9m + C98 * im3_8m -
              C96 * im4_7m - C92 * im5_6m;
    vec[15] = im1 - C94 * re2_9m - C98 * re3_8m + C96 * re4_7m + C92 * re5_6m + C93 * im2_9p - C97 * im3_8p -
              C95 * im4_7p + C91 * im5_6p;
    vec[16] = re1 + C91 * re2_9p + C93 * re3_8p - C95 * re4_7p - C97 * re5_6p + C92 * im2_9m + C94 * im3_8m +
              C96 * im4_7m + C98 * im5_6m;
    vec[17] = im1 - C92 * re2_9m - C94 * re3_8m - C96 * re4_7m - C98 * re5_6m + C91 * im2_9p + C93 * im3_8p -
              C95 * im4_7p - C97 * im5_6p;
}
