/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "peaq.h"
#include "util.h"

#define NUM_BANDS 40                          // filterbank bands
#define SUBSAMP_FB 16                         // filterbank subsampling
#define SUBSAMP_EP 6                          // exitation pattern subsampling
#define SUBSAMP_TOT (SUBSAMP_FB * SUBSAMP_EP) // pcm -> excitation pattern subsampling

// constants used by filters influcenced by subsampling
#define SUBSAMP_TIME SUBSAMP_FB                // can be 32 for quicker fallof
#define SPREADING_C ((double)SUBSAMP_TIME)     // frequency domain spreading
#define SMEARING_TAPS (12 * 32 / SUBSAMP_TIME) // time domain smearing 1 filter legnth
#define SMEARING_OFF (SMEARING_TAPS / 2 - 1)   // time domain smearing 1
#define SMEARING_C (6.0 * SUBSAMP_TIME)        // time domain smearing 2

#define CHUNK_SIZE_EP 20                            // processing chunk
#define CHUNK_SIZE_FB (CHUNK_SIZE_EP * SUBSAMP_EP)  // filterbank samples
#define CHUNK_SIZE_PCM (CHUNK_SIZE_FB * SUBSAMP_FB) // pcm samples

typedef struct {
    double* coef;
    double* taps;
    int     len;
} FirFilter;

typedef struct {
    Complex z1;
    Complex a;
    Complex b;
    Complex w;
} RfirFilter;

typedef struct {
    RfirFilter rfir[3];
    float      zn[CHUNK_SIZE_PCM];
    int        len;
} FdcFilter;

typedef struct {
    FdcFilter fdc;
    float     pcm_buf[CHUNK_SIZE_PCM];
    int       delay;
} PeaqBandFilter;

struct peaq {
    double         scaling_fac;             // playback level
    float          pcm_buf[CHUNK_SIZE_PCM]; // input pcm data
    int            pcm_buf_fill;            // number of samples in pcm_buf
    PeaqBandFilter filter_bank[NUM_BANDS];
    double         re[NUM_BANDS][CHUNK_SIZE_FB]; // real filterbank output
    double         im[NUM_BANDS][CHUNK_SIZE_FB]; // imag filterbank output
    double         upper_spread[NUM_BANDS];      // frequency domain spreading filter
    FirFilter      smearing_1_filter[NUM_BANDS];
    double         smearing_2_filter[NUM_BANDS];
    double         ear_weight[NUM_BANDS];
    double         internal_noise[NUM_BANDS];
    double         forward_mask[NUM_BANDS];
    double         loudness_e[NUM_BANDS]; // excitation at threshold
    double         loudness_s[NUM_BANDS]; // threshold index
    double*        loudness;              // specific loudness buffer, column major
    int            loudness_size;         // size of loudness buffer / NUM_BANDS
    int            loudness_frames;       // number of frames in loudness buffer
};

typedef struct {
    float freq;
    int   length;
    int   delay;
} FilterBand;

static const FilterBand FILTER_BANDS[NUM_BANDS] = {
    {50.00, 1456, 1},    {116.19, 1438, 10},  {183.57, 1406, 26},  {252.82, 1362, 48},  {324.64, 1308, 75},
    {399.79, 1244, 107}, {479.01, 1176, 141}, {563.11, 1104, 177}, {652.97, 1030, 214}, {749.48, 956, 251},
    {853.65, 884, 287},  {966.52, 814, 322},  {1089.25, 748, 355}, {1223.10, 686, 386}, {1369.43, 626, 416},
    {1529.73, 570, 444}, {1705.64, 520, 469}, {1898.95, 472, 493}, {2111.64, 430, 514}, {2345.88, 390, 534},
    {2604.05, 354, 552}, {2888.79, 320, 569}, {3203.01, 290, 584}, {3549.90, 262, 598}, {3933.02, 238, 610},
    {4356.27, 214, 622}, {4823.97, 194, 632}, {5340.88, 176, 641}, {5912.30, 158, 650}, {6544.03, 144, 657},
    {7242.54, 130, 664}, {8014.95, 118, 670}, {8869.13, 106, 676}, {9813.82, 96, 681},  {10858.63, 86, 686},
    {12014.24, 78, 690}, {13292.44, 70, 694}, {14706.26, 64, 697}, {16270.13, 58, 700}, {18000.02, 52, 703}};


static inline void move_float(const float* src, float* dst, int n) { memmove(dst, src, n * sizeof(float)); }

static inline void move_double(const double* src, double* dst, int n) { memmove(dst, src, n * sizeof(double)); }

static inline void zero_float(float* buf, int n) { memset(buf, 0, n * sizeof(float)); }

// vector multiplication out[] = in[] * fac
static inline void vmul(const double* restrict in, double* restrict out, double fac, int len)
{
    for (int i = 0; i < len; i++) {
        out[i] = in[i] * fac;
    }
}

static void alloc_fir(FirFilter* fir, int len)
{
    fir->coef = calloc(sizeof(fir->coef[0]), len);
    fir->taps = calloc(sizeof(fir->taps[0]), len);
    fir->len  = len;
}

static void free_fir(FirFilter* fir)
{
    free(fir->coef);
    free(fir->taps);
    memset(fir, 0, sizeof(*fir));
}

Peaq* peaq_init(double playback_level)
{
    Peaq* pq = calloc(1, sizeof(*pq));

    pq->scaling_fac = exp10(playback_level / 20.0);

    for (int k = 0; k < NUM_BANDS; k++) {
        const FilterBand f   = FILTER_BANDS[k];
        PeaqBandFilter*  pbf = &pq->filter_bank[k];
        pbf->delay           = f.delay;

        // init FDC filter
        assert(ARRAY_SIZE(pbf->fdc.rfir) == 3);
        for (int i = 0; i < 3; i++) {
            const double gamma = M_PI * f.freq / PEAQ_SAMPLERATE;
            const double omega = 2 * (gamma + (i - 1) * M_PI / f.length);
            const double sigma = 4.0 / ((i == 1 ? -2 : 4) * f.length);
            pbf->fdc.rfir[i].a = cneg(cexpi(f.length * omega));
            pbf->fdc.rfir[i].b = cexpi(omega);
            pbf->fdc.rfir[i].w = cmul(cmplx(-sigma, 0), cexpi(-gamma * f.length));
        }
        pbf->fdc.len = f.length;

        // outer & middle ear weights
        const double fk   = f.freq / 1000.0;
        const double omw  = -0.6 * 3.64 * pow(fk, -0.8) + 6.5 * exp(-0.6 * sqr(fk - 3.3)) - exp10(-3) * pow(fk, 3.6);
        pq->ear_weight[k] = exp10(omw / 20);

        // internal noise
        pq->internal_noise[k] = exp10(0.4 * 0.364 * pow(fk, -0.8));

        // time domain smearing 1 filter
        alloc_fir(&pq->smearing_1_filter[k], SMEARING_TAPS);
        for (int i = 0; i < SMEARING_TAPS; i++) {
            pq->smearing_1_filter[k].coef[i] = sqr(cos((M_PI * (i - SMEARING_OFF)) / SMEARING_TAPS));
        }

        // time domain smearing 2
        // 0.008 + 100.0 / f * (0.05 - 0.008) deviates from std!
        pq->forward_mask[k] = exp(-SMEARING_C / ((0.008 + 100.0 / f.freq * (0.05 - 0.008)) * PEAQ_SAMPLERATE));

        // specific loudness
        pq->loudness_e[k] = exp10(0.364 * pow(f.freq / 1000.0, -0.8));
        pq->loudness_s[k] = exp10(0.1 * (-2 - 2.05 * atan(f.freq / 4000.0) - 0.75 * atan(sqr(f.freq / 1600.0))));
    }

    return pq;
}

void peaq_free(Peaq* pq)
{
    if (pq) {
        free(pq->loudness);
        for (int k = 0; k < NUM_BANDS; k++) {
            free_fir(&pq->smearing_1_filter[k]);
        }
        memset(pq, 0, sizeof(*pq));
        free(pq);
    }
}

static void scale_input(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.3
    for (int i = 0; i < CHUNK_SIZE_PCM; i++) {
        pq->pcm_buf[i] *= pq->scaling_fac;
    }
}

// finite impulse respone filter for time domain smearing with 6x subsampled output
static void tds_fir(FirFilter* fir, const double* restrict input, double* restrict output)
{
    // the subsampling has an offset of SMEARING_TAPS / 2. this function itself has a offset of SUBSAMP_EP,
    // the rest of the offset is handled in append_loudness
    assert(SMEARING_TAPS % SUBSAMP_EP == 0);
    assert(SUBSAMP_EP < fir->len);

    double* restrict coeff = fir->coef;
    double* restrict taps  = fir->taps;

    for (int i = 0; i < CHUNK_SIZE_EP; i++) {
        // shift taps to the right and fill with new samples
        move_double(taps, taps + SUBSAMP_EP, fir->len - SUBSAMP_EP);
        for (int j = 0; j < SUBSAMP_EP; j++) {
            taps[j] = input[i * SUBSAMP_EP + SUBSAMP_EP - 1 - j];
        }
        // filter kernel
        float acc = 0;
        for (int t = 0; t < fir->len; t++) {
            acc += coeff[t] * taps[t];
        }
        output[i] = acc;
    }
}

// recursive finite impulse response filter step, zn is delayed input
static inline void rfir_step(RfirFilter* f, float input, float zn)
{
    Complex x = {input + zn * f->a.r + f->z1.r, zn * f->a.i + f->z1.i};
    f->z1     = (Complex){x.r * f->b.r + x.i * -f->b.i, x.i * f->b.r + x.r * f->b.i};
}

// frequency domain convolution filter with output subsampling
static void
subsamp_fdc(FdcFilter* fdc, int subsamp, const float* restrict input, double* restrict re, double* restrict im, int len)
{
    // see Thilo Thiede, "Perceptual Audio Quality Assessment using a Non-Linear Filter Bank" §3.5.2
    // note there are errors in the formulas concerning signs and scaling in that document
    assert(len % subsamp == 0);
    assert(ARRAY_SIZE(fdc->rfir) == 3);

    float* restrict zn = fdc->zn;
    move_float(input, zn + fdc->len, len - fdc->len);

    RfirFilter f0 = fdc->rfir[0];
    RfirFilter f1 = fdc->rfir[1];
    RfirFilter f2 = fdc->rfir[2];

    for (int i = 0; i < len; i++) {
        // output subsampling
        if (i % subsamp == 0) {
            const Complex x = cadd(cadd(cmul(f0.z1, f0.w), cmul(f1.z1, f1.w)), cmul(f2.z1, f2.w));
            re[i / subsamp] = x.r;
            im[i / subsamp] = x.i;
        }
        // filter kernel
        rfir_step(&f0, input[i], zn[i]);
        rfir_step(&f1, input[i], zn[i]);
        rfir_step(&f2, input[i], zn[i]);
    }

    fdc->rfir[0] = f0;
    fdc->rfir[1] = f1;
    fdc->rfir[2] = f2;

    move_float(input + len - fdc->len, zn, fdc->len);
}

static void band_filter(PeaqBandFilter* pbf, const float* input, double* re, double* im, int len)
{
    assert(pbf->delay < len);
    float* restrict pcm_buf = pbf->pcm_buf;

    // copy samples to work buffer accorting to delay
    move_float(input, pcm_buf + pbf->delay, len - pbf->delay);

    // apply filter bank
    subsamp_fdc(&pbf->fdc, SUBSAMP_FB, pcm_buf, re, im, len);

    // copy delayed samples to work buffer for next call
    move_float(input + len - pbf->delay, pcm_buf, pbf->delay);
}

static void filter_bank(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.5
    // output size is 1/16 due to subsampling
    for (int k = 0; k < NUM_BANDS; k++) {
        band_filter(&pq->filter_bank[k], pq->pcm_buf, pq->re[k], pq->im[k], CHUNK_SIZE_PCM);
    }
}

static void ear_weighting(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.6
    for (int k = 0; k < NUM_BANDS; k++) {
        vmul(pq->re[k], pq->re[k], pq->ear_weight[k], CHUNK_SIZE_FB);
        vmul(pq->im[k], pq->im[k], pq->ear_weight[k], CHUNK_SIZE_FB);
    }
}

static void frequency_domain_spreading(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.7
    const double z0   = 7.0 * asinh(50.0 / 650.0);
    const double z39  = 7.0 * asinh(18000.02 / 650.0);
    const double dist = pow(0.1, (z39 - z0) / 780.0);
    const double a    = exp(-SPREADING_C / 4800.0);
    const double b    = 1 - a;

    double* cu = pq->upper_spread;
    double  re[NUM_BANDS];
    double  im[NUM_BANDS];

    for (int i = 0; i < CHUNK_SIZE_FB; i++) {
        for (int k = 0; k < NUM_BANDS; k++) {
            re[k] = pq->re[k][i];
            im[k] = pq->im[k][i];
        }

        for (int k = 0; k < NUM_BANDS; k++) {
            double d1 = pq->re[k][i];
            double d2 = pq->im[k][i];

            const double fc = FILTER_BANDS[k].freq;
            const double l  = 10.0 * log10(sqr(d1) + sqr(d2));
            const double s  = fmax(4.0, 24.0 + 230.0 / fc - 0.2 * l);

            cu[k] = a * pow(dist, s) + b * cu[k];

            for (int j = k + 1; j < NUM_BANDS; j++) {
                d1 *= cu[k];
                d2 *= cu[k];
                re[j] += d1;
                im[j] += d2;
            }
        }

        const double cl = pow(dist, 31.0);
        double       d1 = 0.0;
        double       d2 = 0.0;

        for (int k = NUM_BANDS - 1; k >= 0; k--) {
            d1           = d1 * cl + re[k];
            d2           = d2 * cl + im[k];
            pq->re[k][i] = d1;
            pq->im[k][i] = d2;
        }
    }
}

static void rectification(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.8
    for (int k = 0; k < NUM_BANDS; k++) {
        double* restrict a_re = pq->re[k];
        double* restrict a_im = pq->im[k];
        for (int i = 0; i < CHUNK_SIZE_FB; i++) {
            a_re[i] = sqr(a_re[i]) + sqr(a_im[i]);
        }
    }
}

static void time_domain_smearing_1(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.9
    // output size is 1/6 due to subsampling
    const double c = 2 * 0.9761 / SMEARING_TAPS;
    for (int k = 0; k < NUM_BANDS; k++) {
        // pq->im[k] is used as temporary buffer
        tds_fir(&pq->smearing_1_filter[k], pq->re[k], pq->im[k]);
        vmul(pq->im[k], pq->re[k], c, CHUNK_SIZE_EP);
    }
}

static void internal_noise(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.10
    for (int k = 0; k < NUM_BANDS; k++) {
        const float noise = pq->internal_noise[k];
        for (int i = 0; i < CHUNK_SIZE_EP; i++) {
            pq->re[k][i] += noise;
        }
    }
}

static void time_domain_smearing_2(Peaq* pq)
{
    // ITU-R BS.1387 §2.2.11
    for (int k = 0; k < NUM_BANDS; k++) {
        const float a = pq->forward_mask[k];
        const float b = 1 - a;
        pq->re[k][0]  = a * pq->smearing_2_filter[k] + b * pq->re[k][0];
        for (int i = 1; i < CHUNK_SIZE_EP; i++) {
            pq->re[k][i] = a * pq->re[k][i - 1] + b * pq->re[k][i];
        }
        pq->smearing_2_filter[k] = pq->re[k][CHUNK_SIZE_EP - 1];
    }
}

static void specific_loudness(Peaq* pq)
{
    // ITU-R BS.1387 §3.3
    for (int k = 0; k < NUM_BANDS; k++) {
        const float e = pq->loudness_e[k];
        const float s = pq->loudness_s[k];
        const float c = 1.26539 * 24.0 / NUM_BANDS * pow((e / 10000.0) / s, 0.23);
        for (int i = 0; i < CHUNK_SIZE_EP; i++) {
            pq->re[k][i] = c * fmax(pow(1.0 - s + s * pq->re[k][i] / e, 0.23) - 1.0, 0);
        }
    }
}

// append specific loudness values to loudness buffer
static void append_loudness(Peaq* pq)
{
    const int offset = pq->loudness_size;
    // skip first value to compensate for TDS 1 subsampling offset, not needed if SMEARING_TAPS == 12
    assert(SMEARING_TAPS == 24);
    const int skip = offset == 0;
    // update counters
    pq->loudness_size += CHUNK_SIZE_EP - skip;
    pq->loudness_frames += (pq->pcm_buf_fill + SUBSAMP_TOT - 1) / SUBSAMP_TOT - skip;
    pq->pcm_buf_fill = 0;
    // grow size of output buffer
    pq->loudness = realloc(pq->loudness, pq->loudness_size * NUM_BANDS * sizeof(*pq->loudness));
    // copy specific loudness values
    for (int k = 0; k < NUM_BANDS; k++) {
        for (int i = skip; i < CHUNK_SIZE_EP; i++) {
            pq->loudness[(offset + i - skip) * NUM_BANDS + k] = pq->re[k][i];
        }
    }
}

static void process(Peaq* pq)
{
    scale_input(pq);
    filter_bank(pq); // 16x subsampling
    ear_weighting(pq);
    frequency_domain_spreading(pq);
    rectification(pq);
    time_domain_smearing_1(pq); // 6x subsampling
    internal_noise(pq);
    time_domain_smearing_2(pq);
    specific_loudness(pq);
    append_loudness(pq);
}

// feed new samples into peaq, there is no restriction on len
void peaq_update(Peaq* pq, const float* input, int len)
{
    int processed = 0;
    int remaining = len;

    while (remaining > 0) {
        const int n = imin(CHUNK_SIZE_PCM - pq->pcm_buf_fill, remaining);
        move_float(&input[processed], &pq->pcm_buf[pq->pcm_buf_fill], n);
        pq->pcm_buf_fill += n;

        if (pq->pcm_buf_fill == CHUNK_SIZE_PCM) {
            process(pq);
        }

        processed += n;
        remaining -= n;
    }
}

// call when there is no more pcm to process
void peaq_finish(Peaq* pq)
{
    zero_float(&pq->pcm_buf[pq->pcm_buf_fill], CHUNK_SIZE_PCM - pq->pcm_buf_fill);
    process(pq);
}

static float calc_mld(Peaq* pq1, Peaq* pq2, int segment_size, FILE* outfile)
{
    const int len           = imin(pq1->loudness_frames, pq2->loudness_frames);
    float     max_loud_diff = 0;

    for (int s = 0; s < len; s += segment_size) {
        double segment_loud_diff = 0;

        for (int i = 0; i < imin(segment_size, len - s); i++) {
            double* restrict loud1     = pq1->loudness + (s + i) * NUM_BANDS;
            double* restrict loud2     = pq2->loudness + (s + i) * NUM_BANDS;
            double           loud_diff = 0;

            for (int k = 0; k < NUM_BANDS; k++) {
                loud_diff += fabs(loud1[k] - loud2[k]);
            }
            segment_loud_diff = fmax(loud_diff, segment_loud_diff);
        }

        if (outfile) {
            fprintf(outfile, "%.7f\n", segment_loud_diff);
        }
        max_loud_diff = fmax(segment_loud_diff, max_loud_diff);
    }

    return max_loud_diff;
}

void peaq_print_mld(Peaq* pq1, Peaq* pq2, int segment_size, FILE* segment_outfile)
{
    calc_mld(pq1, pq2, segment_size, segment_outfile);
}

float peaq_get_mld(Peaq* pq1, Peaq* pq2) { return calc_mld(pq1, pq2, CHUNK_SIZE_EP, NULL); }
