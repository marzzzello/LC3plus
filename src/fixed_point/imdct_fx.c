/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void ProcessingIMDCT(
    Word32       y[],       /* i:   spectra data */
    Word16 *     y_e,       /* i:   spectral data exponent */
    const Word16 w[],       /* i:   window coefficients including normalization of sqrt(2/N) and scaled by 2^4 */
    Word16       mem[],     /* i/o: overlap add memory */
    Word16 *     mem_e,     /* i/o: overlap add exponent */
    Word16       x[],       /* o:   time signal out */
    Word16       wLen,      /* i:   window length */
    Word16       N,         /* i:   block size */
    Word16       memLen,    /* i:   overlap add buffer size */
    Word16       frame_dms, /* i:   frame size in ms */
    Word16     concealMethod,     /* i:   concealment method */
    Word16     bfi,               /* i:   bad frame indicator */
    Word16     prev_bfi,          /* i:   previous bad frame indicator */
    Word16     nbLostFramesInRow, /* i: number of consecutive lost frames */
    AplcSetup *plcAd,             /* i: advanced plc struct */
    Word8 *scratchBuffer)
{
    Counter i;
    Word16  o, z, m, s;
    Word16  y_s, mem_s, max_bw;
    Word32  L_tmp;
    Word32 *workBuffer;

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Word16  o, z, m, s;
        Word16  y_s, mem_s, max_bw;
        Word32  L_tmp;
        Counter i;
        Word32 *workBuffer;
        Word16 mem_i_win;
        Word16 w_taper_win;
    };
    Dyn_Mem_In("ProcessingIMDCT", sizeof(struct _dynmem));
#endif


    test(); test(); test();
    IF (sub(bfi, 1) != 0 || concealMethod == 0 || sub(concealMethod, 4) == 0 || sub(concealMethod, 5) == 0)
    {
        workBuffer = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 4 * MAX_LEN bytes */

        /* Init (constant per sample rate) */
        z      = 2 * N - wLen; /* number of leading zeros in window */
        m      = N >> 1;       /* half block size */
        o      = m - z;
        max_bw = 0;

        SWITCH (frame_dms)
        {
        case 25:
            max_bw = MAX_BW >> 2; move16();
            BREAK;
        case 50:
            max_bw = MAX_BW >> 1; move16();
            BREAK;
        case 100:
            max_bw = MAX_BW; move16();
            BREAK;
        }
        
        if (N > max_bw)
            basop_memset(&y[max_bw], 0, (N - max_bw) * sizeof(*y));

        /* Start Processing */
        y_s = getScaleFactor32_0(y, N);
        IF (sub(y_s, 32) < 0)
        {
            FOR (i = 0; i < N; i++)
            {
                y[i] = L_shl(y[i], y_s);
            }
            *y_e = sub(*y_e, y_s);

            dct_IV(y, y_e, N, workBuffer);

            y_s  = getScaleFactor32(y, N);
            y_s  = sub(y_s, 1);
            *y_e = sub(*y_e, y_s + 3); /* mdct window is scaled by pow(2,x) */
            /* N<=20 only happens for 2.5ms frames in NB */
            if (sub(N, 20) <= 0)
            {
                *y_e = add(*y_e, 2);
            }
            else if (sub(N, 120) <= 0)
            {
                *y_e = add(*y_e, 1);
            }
        }
        ELSE
        {
            *y_e = 0;  move16();
        }

        mem_s = getScaleFactor16_0(mem, memLen);
        IF (sub(mem_s, 16) < 0)
        {
            mem_s  = sub(mem_s, 1);
            *mem_e = sub(*mem_e, mem_s);
        }
        ELSE
        {
            *mem_e = *y_e;  move16();
        }

        s = sub(*mem_e, *y_e);

    IF (s > 0)
    {
        y_s  = sub(y_s, s);
        *y_e = add(*y_e, s);
    }
    ELSE
    {
        mem_s  = add(mem_s, s);
        *mem_e = sub(*mem_e, s);
    }

        mem_s = s_max(mem_s, -31);
        y_s   = s_max(y_s, -31);

        UNUSED(prev_bfi);
        UNUSED(nbLostFramesInRow);
        UNUSED(plcAd);

        { /* regular operation */
            FOR (i = 0; i < o; i++)
            {
                L_tmp =
                    L_sub(L_shl(L_deposit_h(mem[i]), mem_s), Mpy_32_16(L_shl(y[m + i + z], y_s), w[4 * m - 1 - i - z]));
                x[i] = round_fx(L_tmp);  move16();
            }
            FOR (i = 0; i < m; i++)
            {
                L_tmp = L_add(L_shl(L_deposit_h(mem[i + o]), mem_s),
                              Mpy_32_16(L_shl(y[2 * m - 1 - i], y_s), w[3 * m - 1 - i]));
                x[i + o] = round_fx(L_tmp);  move16();
            }
        }

        FOR (i = 0; i < m; i++)
        {
            L_tmp            = L_negate(Mpy_32_16(L_shl(y[i], y_s), w[m - 1 - i]));
            x[3 * m - z + i] = round_fx(L_tmp);  move16();
        }

        FOR (i = 0; i < m; i++)
        {
            L_tmp                = L_negate(Mpy_32_16(L_shl(y[i], y_s), w[m + i]));
            x[3 * m - z - 1 - i] = round_fx(L_tmp);  move16();
        }

        basop_memmove(mem, &x[N], memLen * sizeof(Word16));

        *mem_e = *y_e;  move16();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}
/* End Processing */


void Processing_ITDA_WIN_OLA(
    Word32       L_x_tda[], /* i:     X_TDA buffer data   =  "y"  DCT-IV output */
    Word16 *     y_e,       /* i/o:   x_tda  input exponent "y_e"   ,   x output exponent */
    const Word16 w[],       /* i:     window coefficients including normalization of sqrt(2/N) and scaled by 2^4 */
    Word16       mem[],     /* i/o:  overlap add memory */
    Word16 *     mem_e,     /* i/o:  overlap add exponent */
    Word16       x[],       /* o:   time signal out */
    Word16       wLen,      /* i:   window length */
    Word16       N,         /* i:   block size */
    Word16       memLen     /* i:   overlap add buffer size */
    )
{
    /* Declarations */
    Word16  i, o, z, m, s;
    Word16  y_s, mem_s;
    Word32  L_tmp;
    Word32 *L_y;
    Word16 fs_idx, tmp_w, w_factor;
    Word16 factorITDA[5]= { 25905 ,      18318   ,    22435   ,    25905   ,    31727};

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("Processing_ITDA_WIN_OLA", sizeof(struct {
                   Word16  i, o, z, m, s;
                   Word16  y_s, mem_s;
                   Word32  L_tmp;
                   Word32 *L_y;
               }));
#endif


    /* Init (constants  per sample rate) */
    z = 2 * N - wLen; /* number of leading zeros in window */
    m = N >> 1;       /* half block size */
    o = m - z;


    L_y = L_x_tda; /* use same variables naming as in IMDCT for DCT-IV output  signal y */
   
    y_s = getScaleFactor32(L_y, N);        

    y_s = sub(y_s, 1); /*  add 1 bit margin  , y_s is now initial tda upscaling factor */


    *y_e = sub(add(*y_e,1),y_s); /*  handle W scale down by 2^(3) , as mdct synthesis window  was upscaled by  pow(2,x)  x=2 for NB otherwise 3  */

   

    mem_s = getScaleFactor16_0(mem, memLen);

    IF (sub(mem_s, 16) < 0)
    {
        mem_s  = sub(mem_s, 1);      /* one bit margin */
        *mem_e = sub(*mem_e, mem_s); /*adjusted mem exponent due to new scale */
    }
    ELSE
    {
        *mem_e = 0;  move16();
    }

    s = sub(*mem_e, *y_e); /*  */

    IF (s > 0)
    {
        y_s  = sub(y_s, s);     /*  new , reduced upshift of TDA  in window application  loop */
        *y_e = add(*y_e, s);    /*  resulting new exp y_e  for output signal  */
    }
    ELSE
    {
        mem_s  = add(mem_s, s);   /*  s negative or zero,  new , decreased upshift of OLAmem   in loop */
        *mem_e = sub(*mem_e, s);   /*   resulting new exp mem_e  for OLA_mem output signal  */
    }

    fs_idx = mult(N,(Word16)(32768.0/99.0)); /* truncation needed , i.e no rounding can be applied here */  
    w_factor = factorITDA[fs_idx]; move16();
    
    
    
    FOR (i = 0; i < o; i++)
    {   
        tmp_w =  mult_r( w[4 * m - 1 - i - z], w_factor  );
        L_tmp =  L_sub(L_shl_sat(L_deposit_h(mem[i]), mem_s), 
                       Mpy_32_16(L_shl(L_y[m + i + z], y_s), tmp_w ));
        x[i]  = round_fx_sat(L_tmp);  move16();
    }

    FOR (i = 0; i < m; i++)
    {   
        tmp_w    =  mult_r( w[3 * m - 1 - i] , w_factor  );
        L_tmp    = L_add(L_shl_sat(L_deposit_h(mem[i + o]), mem_s),
                         Mpy_32_16(L_shl(L_y[2 * m - 1 - i], y_s),tmp_w ));
        x[i + o] = round_fx_sat(L_tmp);  move16();
    }

    FOR (i = 0; i < m; i++)
    {
        tmp_w            = mult_r( w[m - 1 - i] , w_factor );
        L_tmp            = L_negate(Mpy_32_16(L_shl(L_y[i], y_s), tmp_w));
        x[3 * m - z + i] = round_fx(L_tmp);  move16();
    }

    FOR (i = 0; i < m; i++)
    {
        tmp_w                =  mult_r( w[m + i] , w_factor  );
        L_tmp                = L_negate(Mpy_32_16(L_shl(L_y[i], y_s), tmp_w ));
        x[3 * m - z - 1 - i] = round_fx(L_tmp);  move16();
    }

    FOR (i = 0; i < memLen; i++)
    {
        mem[i] = x[N + i];  move16();
    }
    *mem_e = *y_e;  move16(); /* set OLA mem  exp to  x_Fx exponent*/



#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

