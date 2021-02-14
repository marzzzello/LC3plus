/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef IIS_FFT_H
#define IIS_FFT_H

#include "../structs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct T_IIS_FFT* HANDLE_IIS_FFT;

typedef enum {
    IIS_FFT_NO_ERROR = 0,
    IIS_FFT_INTERNAL_ERROR,     /**< a mystical error appeard */
    IIS_FFT_LENGTH_ERROR,       /**< the requested fft length is not supported */
    IIS_FFT_MEMORY_ERROR        /**< memory allocation failed */
} IIS_FFT_ERROR;

typedef enum {
    IIS_FFT_FWD = -1,           /**< forward transform */
    IIS_FFT_BWD =  1            /**< inverse / backward transform */
} IIS_FFT_DIR;


/*!
 *  \brief          n-point complex FFT
 *
 *                  There are optimized FFTs for lengths 2, 3, 4, 7, 8, 9, 15, 16, 32, 60, 64, 128, 
 *                  240, 256, 384, 480, 512, 768, 1024. Other lengths below 1024 use a stack allocated 
 *                  buffer and offer reasonable speed. Above 1024 a buffer is allocated each time 
 *                  iis_fftf() is called resulting in reduced performance.
 *
 *                  >>>>>> DO NOT USE UNOPTIMIZED LENGTHS IN PRODUCTION CODE! <<<<<<                          
 *
 *  \param[in,out]  vec     pointer to data, interleaved real / imaginary
 *  \param[in]      length  length of fft (number of real/imaginary pairs)
 *
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_iis_fftf(LC3_FLOAT* vec, LC3_INT length);
/*!
 *  \brief          n-point inverse complex FFT
 *
 *                  The output is not normalized. See iis_fftf() for optimized lengths.
 *      
 *                  >>>>>> DO NOT USE UNOPTIMIZED LENGTHS IN PRODUCTION CODE! <<<<<<                          
 *
 *  \param[in,out]  vec     pointer to data, interleaved real / imaginary
 *  \param[in]      length  length of fft (number of real/imaginary pairs)
 *
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_iis_ifftf(LC3_FLOAT* vec, LC3_INT length);

/*!
 *  \brief          allocate and initialize a new real FFT instance.
 *
 *  \param[in,out]  handle      pointer to FFT handle
 *  \param[in]      len         transform length, must be an even number
 *  \param[in]      sign        IIS_FFT_FWD(-1) for forward, IIS_FFT_BWD(1) for backward transform
 *                              BEWARE OF THE SIGNS!
 *
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_IIS_RFFT_Create(HANDLE_IIS_FFT* handle, LC3_INT len, IIS_FFT_DIR sign);

/*!
 *  \brief          allocate and initialize a new complex FFT instance
 *
 *  \param[in,out]  handle      pointer to FFT handle
 *  \param[in]      len         transform length
 *  \param[in]      sign        IIS_FFT_FWD(-1) for forward, IIS_FFT_BWD(1) for backward transform
 *                              BEWARE OF THE SIGNS !!!!!!
 *
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_IIS_CFFT_Create(HANDLE_IIS_FFT* handle, LC3_INT len, IIS_FFT_DIR sign);

/*!
 *  \brief          computes the forward or backward fourier transform of a real signal
 *
 *                  For complex data (in or out) the real part of the Nyquist band (len / 2 + 1) is stored
 *					in the imaginary part of the DC band (0). This allows for the complex data of real to complex
 *					transforms to fit into the same buffer. For this to work length must be even.
 *
 *					Complex to real transforms are normalized (1.0/len). Input and ouput buffers may be identical.
 *
 *  \param[in]      handle      FFT handle
 *  \param[in]      in          pointer to the input array containing real values for the forward transform (FFT)
 *                              or packed complex values (Perm) for the backward transform (IFFT)
 *  \param[out]     out         pointer to the output array containing real values resulted from the backward transform (IFFT)
 *                              or packed complex (perm) values reulted from the forward transform
 *
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_IIS_FFT_Apply_RFFT(HANDLE_IIS_FFT handle, const LC3_FLOAT* in, LC3_FLOAT* out);

/*!
 *  \brief          compute complex backward or forward FFT
 *
 *                  Input and ouput buffers may be identical. Real/imaginary parts may be interleaved.
 *                  The output is not normalized.
 *
 *  \param[in]      handle      FFT handle
 *  \param[in]      in_re       pointer to the input array containing real parts of the signal for the
 *                              forward transform (FFT) or for the backward transform (IFFT)
 *  \param[in]      in_im       pointer to the input array containing imaginary parts of the signal for
 *                              the forward transform (FFT) or for the backward transform (IFFT)
 *  \param[out]     out_re      pointer to the output array containing real values resulted from the
 *                              forward transform (FFT) or from the backward transform (IFFT)
 *  \param[out]     out_im      pointer to the output array containing imaginary values resulted from
 *                              the forward transform (FFT) or from the backward transform (IFFT)
 *
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_IIS_FFT_Apply_CFFT(HANDLE_IIS_FFT handle, const Complex* input, Complex* output);

/*!
 *  \brief          deallocate a FFT instance (complex or real)
 *  \param[in,out]  handle     pointer to FFT handle, set to NULL if call succeeds
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_IIS_xFFT_Destroy(HANDLE_IIS_FFT* handle);

/*!
 *  \brief          deallocate a real FFT instance
 *  \param[in,out]  handle     pointer to FFT handle, set to NULL if call succeeds
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_IIS_RFFT_Destroy(HANDLE_IIS_FFT* handle);

/*!
 *  \brief          deallocate a complex FFT instance
 *  \param[in,out]  handle     pointer to FFT handle, set to NULL if call succeeds
 *  \return         IIS_FFT_NO_ERROR on success
 */
IIS_FFT_ERROR LC3_IIS_CFFT_Destroy(HANDLE_IIS_FFT* handle);


#ifdef __cplusplus
}
#endif

#endif
