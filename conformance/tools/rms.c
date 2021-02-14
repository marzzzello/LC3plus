/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "tinywavein_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* Global defines */
/* K = 16 bit */
#define RMS_MAX_BUF 1024
#define SCALE_16 (1 << 15)
#define SCALE_24 (1 << 23)
#define MAX_DIFF 0.000061035
#define MAX_RMS -101.1008
#define SEGMENT_LENGTH 320
#define SSNR_LOW_THR -50.0
#define SSNR_HIGH_THR -15.0

/* Function declarations */

static void printResult(char *inputFilename1, char *inputFilename2, int totalSamples1, float diffMax, float rms, float ssnr, int segmentLength, int differentSamples, float maxDiffThr, float maxRmsThr);
static void calculateRms(char *inputFilename1, char *inputFilename2, int *differentSamples, int *totalSamples1, int *totalSamples2, double *rmsOut, float *maxDiffOut);
static void calculateSegmentalSnr(char *inputFilename1, char *inputFilename2, float *ssnrOut);
static int checkRmsReached(float rms);
static void printUsage(void);
static void calculateThr(int k, float *maxDiffThr, float *maxRmsThr);

int main(int ac, char *av[])
{
    char     *inputFilename1 = NULL, *inputFilename2 = NULL;
    int       totalSamples1 = 0, totalSamples2 = 0, differentSamples = 0, k = 0;
    float diffMax = 0, ssnr = 0, maxDiffThr = 0, maxRmsThr = 0;
    double rms = 0;

    if(ac < 3)
    {
        printf("    Not enough input arguments!\n");
        printUsage();
    }
    
    if(ac == 4)
    {
        k = atoi(av[3]);
        if(k <= 0 || k > 16)
        {
            printf("    Parameter k has to be in range of 1 ... 16!\n");
            exit(1);
        } else {
            printf("    Using k = %d as threshold!\n", k);
        }
    }
    
    if (ac > 4)
    {
        printf("    Too many input arguments!\n\n");
        printUsage();
    }

    inputFilename1 = av[1]; /* Reference */
    inputFilename2 = av[2]; /* Codec under test */
    
    if(k != 0)
    {
        calculateThr(k, &maxDiffThr, &maxRmsThr);
    } else {
        maxDiffThr = MAX_DIFF;
        maxRmsThr = MAX_RMS;
    }

    calculateRms(inputFilename1, inputFilename2, &differentSamples, &totalSamples1, &totalSamples2, &rms, &diffMax);

    calculateSegmentalSnr(inputFilename1, inputFilename2, &ssnr);

    printResult(inputFilename1, inputFilename2, totalSamples1, diffMax, rms, ssnr, SEGMENT_LENGTH, differentSamples, maxDiffThr, maxRmsThr);

    return 0;
}

void calculateThr(int k, float *maxDiffThr, float *maxRmsThr)
{
    *maxDiffThr = 1 / (pow(2, (k - 2)));
    
    *maxRmsThr = 20 * log10(pow(2, -(k - 1)) / sqrt(12));
}

void printResult(char *inputFilename1, char *inputFilename2, int totalSamples1, float diffMax, float rms, float ssnr, int segmentLength, int differentSamples, float maxDiffThr, float maxRmsThr)
{
    char *maxDiffReached, *maxRmsReached;
    char *rmsFail, *diffFail;
    int  rmsReached = 0;

    if(diffMax < maxDiffThr)
    {
        maxDiffReached = "REACHED";
        diffFail = "PASSED";
    } else {
        maxDiffReached = "NOT REACHED";
        diffFail = "FAILED";
    }

    if(rms < maxRmsThr)
    {
        maxRmsReached = "REACHED";
        rmsFail = "PASSED";
    } else {
        maxRmsReached = "NOT REACHED";
        rmsFail = "FAILED";
    }

    printf("Comparing files: %s and %s \n\n", inputFilename1, inputFilename2);

    printf("    Number of samples compared          : %d\n", totalSamples1);
    printf("    Number of different samples         : %d\n", differentSamples);

    if(differentSamples == 0)
    {
        printf("\n    Input files match exactly!\n\n");
    }

    if(differentSamples != 0)
    {
        printf("    Maximum difference                  : %e ---- %s ---- (threshold for 16-bit resolution is %e)\n", diffMax, maxDiffReached, maxDiffThr);
        printf("    Overall RMS value                   : %f dB ---- %s ---- (threshold for 16-bit resolution is %f dB)\n", rms, maxRmsReached, maxRmsThr);
        printf("    Average SNR value                   : %f dB (%d samples per segment)\n\n", ssnr, segmentLength);

        printf("---- Test on RMS criteria               : %s\n", rmsFail);
        printf("---- Test on max. abs. diff criteria    : %s\n", diffFail);

        rmsReached = checkRmsReached(rms);
        printf("---- Reached RMS criteria               : %d bit\n", rmsReached);
    }

}

void calculateRms(char *inputFilename1, char *inputFilename2, int *differentSamples, int *totalSamples1, int *totalSamples2, double *rmsOut, float *maxDiffOut)
{
    /* Calculate RMS */

    int nSamples1 = 0, nSamples2 = 0, scale = 0, i;
    int sampleRate1 = 0, sampleRate2 = 0, nLength1 = 0, nLength2 = 0, nSamplesRead1 = 0, nSamplesRead2 = 0;
    float sample_buf1_scaled[RMS_MAX_BUF], sample_buf2_scaled[RMS_MAX_BUF];
    int sample_buf1[RMS_MAX_BUF], sample_buf2[RMS_MAX_BUF];
    WAVEFILEIN *in_file1, *in_file2;
    int nChannels1 = 0, nChannels2 = 0, bipsIn1 = 0, bipsIn2 = 0;
    float diffMax = 0.0, rms = 0.0;

    in_file1 = OpenWav(inputFilename1, &sampleRate1, &nChannels1, &nLength1, &bipsIn1);
    in_file2 = OpenWav(inputFilename2, &sampleRate2, &nChannels2, &nLength2, &bipsIn2);

    if(in_file1 == NULL || in_file2 == NULL)
    {
        printf("Error opening wave files!\n");
        exit(1);
    }

    if(bipsIn1 == 16)
    {
        scale = SCALE_16;
    } else if (bipsIn1 == 24)
    {
        scale = SCALE_24;
    } else {
        printf("Bits per sample of input files is not supported!\n");
        exit(1);
    }

    nSamples1 = sampleRate1 / 100;
    nSamples2 = sampleRate2 / 100;

    while(1)
    {
        ReadWavShort(in_file1, sample_buf1, nSamples1 * nChannels1, &nSamplesRead1);
        ReadWavShort(in_file2, sample_buf2, nSamples2 * nChannels2, &nSamplesRead2);

        for(i = 0; i < (int) nSamplesRead1; i++)
        {
            sample_buf1_scaled[i] = (float)sample_buf1[i] / scale;
            sample_buf2_scaled[i] = (float)sample_buf2[i] / scale;
        }

        for(i = 0; i < nSamplesRead1; i++)
        {
            /* Get maximum difference */
            if(fabsf((sample_buf1_scaled[i] - sample_buf2_scaled[i])) > diffMax)
            {
                diffMax = fabsf((sample_buf1_scaled[i] - sample_buf2_scaled[i]));
            }

            if((sample_buf1_scaled[i] - sample_buf2_scaled[i]) != 0)
            {
                rms += (sample_buf1_scaled[i] - sample_buf2_scaled[i]) * (sample_buf1_scaled[i] - sample_buf2_scaled[i]);
                *differentSamples = *differentSamples + 1;
            }
        }

        if (nSamplesRead1 != (nSamples1 * nChannels1))
        {
            *totalSamples1 = *totalSamples1 + nSamplesRead1 * nChannels1;
            *totalSamples2 = *totalSamples2 + nSamplesRead2 * nChannels2;
            break;
        }

        *totalSamples1 = *totalSamples1 + nSamples1 * nChannels1;
        *totalSamples2 = *totalSamples2 + nSamples2 * nChannels2;
    }

    rms = rms / *totalSamples1;
    rms = sqrt(rms);
    rms = 20.0 * log10(rms);

    *rmsOut = rms;
    *maxDiffOut = diffMax;

    CloseWavIn(in_file1);
    CloseWavIn(in_file2);
}

void calculateSegmentalSnr(char *inputFilename1, char *inputFilename2, float *ssnrOut)
{
    float nom = 0, denom = 0, pow1 = 0, ss = 0, ssnr = 0;
    int skip = 0;
    float sample_buf1_scaled[RMS_MAX_BUF], sample_buf2_scaled[RMS_MAX_BUF];
    int sample_buf1[RMS_MAX_BUF], sample_buf2[RMS_MAX_BUF];
    WAVEFILEIN *in_file1, *in_file2;
    int nChannels1 = 0, nChannels2 = 0, bipsIn1 = 0, bipsIn2 = 0;
    int scale = 0, segmentLength = 0, nSegments = 0, i;
    int sampleRate1 = 0, sampleRate2 = 0, nLength1 = 0, nLength2 = 0, nSamplesRead1 = 0, nSamplesRead2 = 0;
    in_file1 = OpenWav(inputFilename1, &sampleRate1, &nChannels1, &nLength1, &bipsIn1);
    in_file2 = OpenWav(inputFilename2, &sampleRate2, &nChannels2, &nLength2, &bipsIn2);

    if(in_file1 == NULL || in_file2 == NULL)
    {
        printf("Error opening wave files!\n");
        exit(1);
    }

    if(bipsIn1 == 16)
    {
        scale = SCALE_16;
    } else if (bipsIn1 == 24)
    {
        scale = SCALE_24;
    } else {
        printf("Bits per sample of input files is not supported!\n");
        exit(1);
    }

    segmentLength = SEGMENT_LENGTH;

    while(1)
    {
        nom = 0;
        denom = 0;
        skip = 0;
        pow1 = 0;
        ReadWavShort(in_file1, sample_buf1, segmentLength, &nSamplesRead1);
        ReadWavShort(in_file2, sample_buf2, segmentLength, &nSamplesRead2);

        for(i = 0; i < (int) nSamplesRead1; i++)
        {
            sample_buf1_scaled[i] = (float)sample_buf1[i] / scale;
            sample_buf2_scaled[i] = (float)sample_buf2[i] / scale;
        }

        /* Check if singal power is in the range of [-50 dB ... -15 dB] */
        for(i = 0; i < nSamplesRead1; i++)
        {
            pow1 += (sample_buf1_scaled[i] * sample_buf1_scaled[i]);
        }

        pow1 = 10 * log10(pow1/ (float) nSamplesRead1);

        if(pow1 < SSNR_LOW_THR || pow1 > SSNR_HIGH_THR)
        {
            skip = 1;
        }

        if(skip == 0)
        {
            nSegments++;
            for(i = 0; i < nSamplesRead1; i++)
            {
                  nom += sample_buf1_scaled[i] * sample_buf1_scaled[i];
                  denom += (sample_buf1_scaled[i] - sample_buf2_scaled[i]) * (sample_buf1_scaled[i] - sample_buf2_scaled[i]);
            }

            denom += (float) nSamplesRead1 * (float) pow(10, -13);

            ss = log10(1.0 + nom/denom);
            ssnr += ss;
        }

        if (nSamplesRead1 != (segmentLength * nChannels1))
        {
            break;
        }
    }

    ssnr = ssnr / (float) nSegments;
    ssnr = 10.0 * log10(pow(10, ssnr) - 1.0);

    *ssnrOut = ssnr;

    CloseWavIn(in_file1);
    CloseWavIn(in_file2);
}

int checkRmsReached(float rms)
{
    int i;
    float currentThr;

    for(i = 16; i > 0; i--)
    {
        currentThr = 20.0 * log10(pow(2, -(i - 1)) / sqrt(12.0));
        if(rms < currentThr)
        {
            break;
        }
    }

    return i;
}

void printUsage(void)
{
    printf("    RMS tool to calculate RMS, max. abs. difference and segmental SNR value between to wave files.\n");
    printf("    Usage: rms ref.wav test.wav [k]\n");
    printf("    The test is done in regards to k = 16, i.e. a 16-bit resolution as default. The user can set k to a range from 1 ... 16 bits. The segment length is set 320 samples.\n");
    
    exit(0);
}

