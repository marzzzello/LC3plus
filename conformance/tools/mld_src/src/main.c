/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "peaq.h"
#include "wav_io.h"
#include "util.h"

#define VERSION "0.3.1"
#define BUF_SIZE 1920

// clang-format off
static const char* HELP_MESSAGE =
    // help lines must not be longer than this------------------------------------>|
    "Maximum Loudness Difference Tool "VERSION"\n"
    "Copyright Fraunhofer IIS 2019\n"
    "Syntax: mld [options] <input1> <input2>\n"
    "options:\n"
    "    -d       Disable multithreading\n"
    "    -h       Print help\n"
    "    -l LEVEL Playback level in dBspl (default: 92)\n"
    "    -o FILE  Write output to file\n"
    "    -s       Print segment values\n"
    "    -v       Print version";

static const char* WAV_ERROR_MESSAGE[] = {
    [WAV_OPEN_ERROR]   = "Failed to open file",
    [WAV_READ_ERROR]   = "File read error",
    [WAV_WRITE_ERROR]  = "File write error",
    [WAV_SYNTAX_ERROR] = "File not a wav",
    [WAV_FORMAT_ERROR] = "Unsupported wav format",
    [WAV_NULL_ERROR]   = "Function called with NULL"
};
// clang-format on

typedef struct {
    bool   print_help;
    bool   print_segments;
    bool   print_version;
    bool   single_thread;
    double playback_level;
    char*  infile1;
    char*  infile2;
    char*  outfile;
} Arguments;

typedef struct {
    Peaq* peaq;
    Wav   wav;
} Process;

// check condition. if it fails, print message to stderr and exit
static void exit_if(bool condition, const char* message)
{
    if (condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static const char* parse_args(Arguments* args, int argc, char** argv)
{
    *args = (Arguments){0};

    for (int i = 1, n = 1; n < argc;) {
        int ninc = 1;
        if (argv[n][0] == '-') {
            // parse flags
            char  flag = argv[n][i++];
            char* parm = argv[n][i] ? &argv[n][i] : argv[n + 1];
            if (strchr("hdsv", flag)) { // flags with no parameter
                ninc = !argv[n][i];
            } else if (strchr("lo", flag)) { // flags with parameter
                ninc += !argv[n][i];
            } else {
                return "Invalid option";
            }
            // set argument values
            args->print_help |= flag == 'h';
            args->print_segments |= flag == 's';
            args->print_version |= flag == 'v';
            args->single_thread |= flag == 'd';
            if (flag == 'l' && !atof_c(parm, &args->playback_level)) {
                return "Expecting float argument";
            } else if (flag == 'o') {
                args->outfile = parm;
            }
        } else {
            // positional arguments
            if (!args->infile1) {
                args->infile1 = argv[n];
            } else if (!args->infile2) {
                args->infile2 = argv[n];
            } else {
                return "Unexpected argument";
            }
        }
        i = ninc ? 1 : i;
        n += ninc;
    }

    if (!args->infile1 || !args->infile2) {
        return "Expecting file argument";
    }
    return NULL;
}

// process a complete wav
static void process(Peaq* peaq, Wav* wav)
{
    float  pcm_data[BUF_SIZE];
    float* pcm_buf[1] = {pcm_data};

    int processed = 0;
    while (processed < wav->length) {
        int nsamp = wav_read(wav, pcm_buf, BUF_SIZE);
        if (nsamp <= 0)
            break;
        peaq_update(peaq, pcm_buf[0], nsamp);
        processed += nsamp;
    }

    peaq_finish(peaq);
}

// thread function for process
static void* process_thread(void* arg)
{
    Process* p = arg;
    process(p->peaq, &p->wav);
    return NULL;
}

int main(int argc, char** argv)
{
    // parse commandline
    Arguments   args   = {0};
    const char* argerr = parse_args(&args, argc, argv);
    if (args.print_help || args.print_version) {
        puts(args.print_help ? HELP_MESSAGE : VERSION);
        exit(EXIT_SUCCESS);
    }
    exit_if(argerr, argerr);

    Process pr1 = {0};
    Process pr2 = {0};
    int     err = 0;

    // open input wavs
    err = wav_reader(&pr1.wav, args.infile1);
    exit_if(err, WAV_ERROR_MESSAGE[err]);
    exit_if(pr1.wav.channels != 1, "Input must be mono");
    exit_if(pr1.wav.samplerate != PEAQ_SAMPLERATE, "Input must be 48 kHz");

    err = wav_reader(&pr2.wav, args.infile2);
    exit_if(err, WAV_ERROR_MESSAGE[err]);
    exit_if(pr2.wav.channels != 1, "Input must be mono");
    exit_if(pr2.wav.samplerate != PEAQ_SAMPLERATE, "Input must be 48 kHz");

    // initialize peaq
    const double playback_level = args.playback_level ? args.playback_level : 92;
    pr1.peaq                    = peaq_init(playback_level);
    pr2.peaq                    = peaq_init(playback_level);

    // process input
    if (args.single_thread) {
        process(pr1.peaq, &pr1.wav);
        process(pr2.peaq, &pr2.wav);
    } else {
        pthread_t thread1, thread2;
        err |= pthread_create(&thread1, NULL, process_thread, &pr1);
        err |= pthread_create(&thread2, NULL, process_thread, &pr2);
        exit_if(err, "Failed to create thread");
        err |= pthread_join(thread1, NULL);
        err |= pthread_join(thread2, NULL);
        exit_if(err, "Failed to join thread");
    }

    // open output file if requested
    FILE* fout = stdout;
    if (args.outfile) {
        fout = fopen(args.outfile, "w");
        exit_if(fout == NULL, "Failed to create output file");
    }

    // print result
    if (args.print_segments) {
        peaq_print_mld(pr1.peaq, pr2.peaq, 10, fout);
    } else {
        fprintf(fout, "maximum loudness difference: %f\n", peaq_get_mld(pr1.peaq, pr2.peaq));
    }

    wav_close(&pr1.wav);
    wav_close(&pr2.wav);
    peaq_free(pr1.peaq);
    peaq_free(pr2.peaq);
    fclose(fout);
}
