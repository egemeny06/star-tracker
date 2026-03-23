#include "simulation.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/* Simple xorshift32 PRNG - embedded friendly, no stdlib rand */
static uint32_t xorshift32(uint32_t *state)
{
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

/* Box-Muller transform: generate Gaussian from uniform */
static float rand_gaussian(uint32_t *state)
{
    float u1 = (float)(xorshift32(state) & 0x7FFFFFFF) / (float)0x7FFFFFFF;
    float u2 = (float)(xorshift32(state) & 0x7FFFFFFF) / (float)0x7FFFFFFF;

    if (u1 < 1e-10f) u1 = 1e-10f;

    return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
}

void add_centroid_noise(CentroidList *list, float noise_sigma, uint32_t seed)
{
    uint32_t state = seed;
    if (state == 0) state = 12345;

    for (uint16_t i = 0; i < list->count; i++) {
        float nx = rand_gaussian(&state) * noise_sigma;
        float ny = rand_gaussian(&state) * noise_sigma;

        list->centroids[i].x += nx;
        list->centroids[i].y += ny;
        list->centroids[i].x_err = noise_sigma;
        list->centroids[i].y_err = noise_sigma;
    }
}

float compute_centroid_rms_error(const CentroidList *noisy,
                                const CentroidList *ground_truth)
{
    uint16_t n = noisy->count;
    if (n > ground_truth->count) n = ground_truth->count;
    if (n == 0) return 0.0f;

    float sum_sq = 0.0f;
    for (uint16_t i = 0; i < n; i++) {
        float dx = noisy->centroids[i].x - ground_truth->centroids[i].x;
        float dy = noisy->centroids[i].y - ground_truth->centroids[i].y;
        sum_sq += dx * dx + dy * dy;
    }

    return sqrtf(sum_sq / (float)n);
}

void compute_centroid_errors(const CentroidList *noisy,
                             const CentroidList *ground_truth,
                             float *errors_x, float *errors_y,
                             uint16_t *count)
{
    uint16_t n = noisy->count;
    if (n > ground_truth->count) n = ground_truth->count;
    *count = n;

    for (uint16_t i = 0; i < n; i++) {
        errors_x[i] = noisy->centroids[i].x - ground_truth->centroids[i].x;
        errors_y[i] = noisy->centroids[i].y - ground_truth->centroids[i].y;
    }
}

void simulation_result_print(const SimulationResult *r)
{
    printf("=== Simulation Result ===\n");
    printf("  Noise sigma:        %.4f pixels\n", r->noise_sigma);
    printf("  Centroid RMS error: %.4f pixels\n", r->centroid_rms);
    printf("  Stars identified:   %d\n", r->stars_identified);
    printf("  Attitude error:     %.6f deg (%.2f arcsec)\n",
           r->attitude_error_deg, r->attitude_error_arcsec);
    printf("  Roll  error:        %.6f deg\n", r->roll_error);
    printf("  Pitch error:        %.6f deg\n", r->pitch_error);
    printf("  Yaw   error:        %.6f deg\n", r->yaw_error);
    printf("=========================\n");
}
