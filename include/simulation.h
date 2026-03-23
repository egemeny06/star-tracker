#ifndef SIMULATION_H
#define SIMULATION_H

#include "centroid.h"
#include "attitude.h"

/*
 * Add Gaussian noise to centroid positions.
 * noise_sigma: standard deviation in pixels (sub-pixel, e.g. 0.1)
 * seed: random seed for reproducibility
 */
void add_centroid_noise(CentroidList *list, float noise_sigma, uint32_t seed);

/*
 * Compute centroid position error statistics between noisy and ground truth.
 * Returns RMS error in pixels.
 */
float compute_centroid_rms_error(const CentroidList *noisy,
                                const CentroidList *ground_truth);

/*
 * Compute per-centroid position error.
 * errors_out must have at least min(noisy->count, truth->count) elements.
 */
void compute_centroid_errors(const CentroidList *noisy,
                             const CentroidList *ground_truth,
                             float *errors_x, float *errors_y,
                             uint16_t *count);

/*
 * Full simulation pipeline:
 * - Given ground truth centroids and noise sigma
 * - Add noise, run star ID, compute attitude
 * - Compare with true attitude
 * - Report attitude error in arcseconds
 */
typedef struct {
    float noise_sigma;          /* input noise (pixels) */
    float centroid_rms;         /* centroid RMS error (pixels) */
    float attitude_error_deg;   /* attitude error (degrees) */
    float attitude_error_arcsec;/* attitude error (arcseconds) */
    float roll_error;           /* roll error (degrees) */
    float pitch_error;          /* pitch error (degrees) */
    float yaw_error;            /* yaw error (degrees) */
    int stars_identified;       /* number of stars matched */
} SimulationResult;

void simulation_result_print(const SimulationResult *r);

#endif /* SIMULATION_H */
