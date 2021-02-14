/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define G192_GOOD_FRAME 0x6B21
#define G192_BAD_FRAME 0x6B20
#define G192_REDUNDANCY_FRAME 0x6B22
#define G192_ZERO 0x007F
#define G192_ONE 0x0081
#define MAX_FLIPS 1024

static FILE *output_bitstream;
static FILE *input_bitstream;
static uint64_t rng_state;

static int read_bitstream_frame_G192(FILE *bitstream_file, int size, int16_t *bytes);
static void write_bitstream_frame_G192(FILE *bitstream_file, int16_t *bits, int size);
static void printUsage();
static void flip(int16_t *bits, int length, int flips);
static double my_rand();

/* Usage: flipG192 bs_in bs_out flips frames seed
 * with bs_in     : G192 input bitstream
 *      bs_out    : G192 output bitstream
 *      flips     : number of flips per frames
 *      frames    : number of frames to be flipped in percent, e.g. 50 for 50 %
 *      seed      : seed for random generator
 *      verbose   : verbose output
 */

int main(int ac, char **av)
{
    char  *inputFilename;
    char  *outputFilename;
    int16_t  bits[2000 * 8];
    int ret = 0, flips = 0, seedLoc = 0, frames = 0, randFrame = 0, frame_counter = 1;
    int totalFr = 0, totalFl = 0, totalFrames = 0, verbose = 0;

    if (ac != 7)
    {
        printUsage();
    }
    
    inputFilename = av[1];
    outputFilename = av[2];
    flips = atoi(av[3]);
    frames = atoi(av[4]) * 100;
    seedLoc = atoi(av[5]);
    verbose = atoi(av[6]);
    
    rng_state = seedLoc;
    
    assert(flips >= 0 && "Number of flips must be >= 0!");
    assert(frames >= 0 && "Number of frames must be >= 0!");
    assert((verbose == 0 || verbose == 1) && "Verbose argument either 0 or 1!");
    assert((inputFilename[0] != '\0') && "Input file name must not be empty!");
    assert((outputFilename[0] != '\0') && "Output file name must not be empty!");
    
    input_bitstream = fopen(inputFilename, "rb");
    output_bitstream = fopen(outputFilename, "wb");
    
    assert(input_bitstream && "Could not open input file!");
    assert(output_bitstream && "Could not open output file!");
    
    while (1)
    {
        ret = read_bitstream_frame_G192(input_bitstream, sizeof(bits), bits);
        if (ret < 0)
        {
            break;
        }
        
        totalFrames++;
        
        randFrame = (int) round(my_rand() * 10000);
        randFrame = frames >= randFrame;
        
        if (flips && randFrame)
        {
            totalFr++;
            if (verbose)
            {
                printf("Braking frame %d ...\n", frame_counter);
            }
            
            flip(bits, ret, flips);
            totalFl += flips;
        }
        
        write_bitstream_frame_G192(output_bitstream, bits, ret);
        frame_counter++;
    }
    
    if (verbose)
    {
        printf("\nFlipped %d frames out of %d and totally %d bits\n", totalFr, totalFrames, totalFl);
    }
}

static void flip(int16_t *bits, int length, int flips)
{
    int i = 0, j = 0;
    int allFlips[MAX_FLIPS] = {0};
    assert(flips <= MAX_FLIPS);

    /* Create array containing bit positions that we need to flip */
    while(i < flips){
        int toFlip = (int) round(my_rand() * length);

        for (j = 0; j < i; j++)
        {
            if (allFlips[j] == toFlip)
            {
                break;
            }
        }
        
        if (j == i)
        {
            allFlips[i++] = toFlip;
        }
    }
    
    for (i = 0; i < flips; i++)
    {
        int value = bits[allFlips[i]];
        
        if (value == G192_ZERO)
        {
            bits[allFlips[i]] = G192_ONE;
        } else {
            bits[allFlips[i]] = G192_ZERO;
        }
    }
}

static int read_bitstream_frame_G192(FILE *bitstream_file, int size, int16_t *bits)
{
    int      i = 0, read = 0;
    uint16_t nbits      = 0;
    int16_t  currentBit = 0, frameIndicator = 0;

    /* Read frame indicator info -> good/bad/redundancy frame */
    read = (int)fread(&frameIndicator, sizeof(frameIndicator), 1, bitstream_file);
    
    if (read != 1)
    {
        return -1;
    }

    /* Read length info */
    read = (int)fread(&nbits, sizeof(nbits), 1, bitstream_file);

    assert(((frameIndicator == G192_GOOD_FRAME) || (frameIndicator == G192_BAD_FRAME) || (frameIndicator == G192_REDUNDANCY_FRAME)) &&
            "Wrong G192 format detected in bitstream file! The sync word could not be recognized!");

    for (i = 0; i < nbits; i++)
    {
        read = (int)fread(&currentBit, sizeof(currentBit), 1, bitstream_file);
        bits[i] = currentBit;
    }

    return nbits;
}

static void write_bitstream_frame_G192(FILE *bitstream_file, int16_t *bits, int size)
{
    int      i           = 0;
    int16_t   currentBit = 0, syncWord = 0;

    /* Write good/bad frame info -> encoder writes only good frames */
    syncWord = G192_GOOD_FRAME;
    fwrite(&syncWord, sizeof(int16_t), 1, bitstream_file);

    /* Write length info */
    fwrite(&size, sizeof(uint16_t), 1, bitstream_file);

    for (i = 0; i < size; i++)
    {
        currentBit = bits[i];
        fwrite(&currentBit, sizeof(currentBit), 1, bitstream_file);
    }
}

static void printUsage()
{
    printf("Usage: flipG192 bs_in bs_out FLIPS FRAMES SEED verbose\n");
    printf("       bs_in    : input bitstream in G192 format\n");
    printf("       bs_out   : input bitstream in G192 format\n");
    printf("       FLIPS    : number of bits to flip in a frame\n");
    printf("       FRAMES   : number of frames to flip in percent, e.g. 50 for 50%%\n");
    printf("       SEED     : seed for random generator\n");
    printf("       verbose  : (0): do not print detailed information, (1) print detailed information \n");
    
    exit(0);
}

static double my_rand()
{
	uint64_t z = (rng_state += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	z = z ^ (z >> 31);
    
    return (z >> 11) * (1.0 / ((uint64_t) 1 << 53));
}

