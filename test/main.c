#include <stdio.h>
#include <math.h>
#include <string.h>
#include "centroid.h"
#include "catalog.h"
#include "star_tracker.h"
#include "attitude.h"
#include "projection.h"
#include "star_id.h"
#include "simulation.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/* Build a known rotation matrix from Euler angles (deg) */
static void euler_to_dcm(float roll_deg, float pitch_deg, float yaw_deg, DCM *dcm)
{
    float r = roll_deg  * (float)(M_PI / 180.0);
    float p = pitch_deg * (float)(M_PI / 180.0);
    float y = yaw_deg   * (float)(M_PI / 180.0);

    float cr = cosf(r), sr = sinf(r);
    float cp = cosf(p), sp = sinf(p);
    float cy = cosf(y), sy = sinf(y);

    dcm->m[0][0] = cy * cp;
    dcm->m[0][1] = cy * sp * sr - sy * cr;
    dcm->m[0][2] = cy * sp * cr + sy * sr;
    dcm->m[1][0] = sy * cp;
    dcm->m[1][1] = sy * sp * sr + cy * cr;
    dcm->m[1][2] = sy * sp * cr - cy * sr;
    dcm->m[2][0] = -sp;
    dcm->m[2][1] = cp * sr;
    dcm->m[2][2] = cp * cr;
}

/* Static allocations for large structures */
static StarCatalog catalog;
static StarPairDB pair_db;

int main(void)
{
    printf("****************************************************\n");
    printf("*  STAR TRACKER - FULL PIPELINE SIMULATION          *\n");
    printf("****************************************************\n\n");

    /* ========================================================= */
    /*  STEP 1: Camera setup                                      */
    /* ========================================================= */
    CameraParams cam;
    camera_init(&cam, 16.0f, 0.0055f, 1024, 1024);
    printf("[1] CAMERA\n");
    printf("    Focal: %.1f mm | Pixel: %.4f mm | %ux%u\n",
           cam.focal_length, cam.pixel_size, cam.width, cam.height);
    printf("    FOV: %.2f x %.2f deg\n\n", cam.fov_x, cam.fov_y);

    /* ========================================================= */
    /*  STEP 2: Load Hipparcos catalog                            */
    /* ========================================================= */
    printf("[2] HIPPARCOS CATALOG\n");
    int loaded = catalog_load_csv(&catalog, "data/hipparcos_catalog.csv");
    if (loaded <= 0) {
        printf("    ERROR: Cannot load catalog!\n");
        return 1;
    }
    printf("    Loaded: %d stars\n", loaded);

    CatalogStar bright[20];
    uint16_t nbright = catalog_filter_by_magnitude(&catalog, 2.0f, bright, 20);
    printf("    Stars brighter than mag 2.0: %u\n\n", nbright);

    /* ========================================================= */
    /*  STEP 3: Define TRUE attitude (boresight pointing)         */
    /* ========================================================= */
    /* Boresight pointing toward Orion (RA~85°, Dec~0°) */
    float true_roll  = 5.0f;
    float true_pitch = 0.0f;
    float true_yaw   = 85.0f;

    DCM true_dcm;
    euler_to_dcm(true_roll, true_pitch, true_yaw, &true_dcm);

    Quaternion true_q;
    dcm_to_quaternion(&true_dcm, &true_q);

    printf("[3] TRUE ATTITUDE\n");
    printf("    Roll: %+.2f | Pitch: %+.2f | Yaw: %+.2f deg\n", true_roll, true_pitch, true_yaw);
    printf("    Quaternion: [%.6f, %.6f, %.6f, %.6f]\n\n",
           true_q.q1, true_q.q2, true_q.q3, true_q.q4);

    /* ========================================================= */
    /*  STEP 4: Project stars to image → ground truth centroids   */
    /* ========================================================= */
    CentroidList ground_truth;
    int nstars = project_stars_to_image(&cam, &true_dcm, &catalog, 6.0f, &ground_truth);

    printf("[4] STAR PROJECTION (Ground Truth)\n");
    printf("    Stars in FOV (mag <= 6.0): %d\n", nstars);
    if (nstars > 0) {
        centroid_list_sort_by_flux(&ground_truth);
        printf("    Brightest centroids (sub-pixel ground truth):\n");
        printf("    %-4s  %-12s  %-12s  %-12s\n", "ID", "X (px)", "Y (px)", "Flux");
        int show = nstars > 8 ? 8 : nstars;
        for (int i = 0; i < show; i++) {
            printf("    %-4u  %-12.4f  %-12.4f  %-12.1f\n",
                   ground_truth.centroids[i].id,
                   ground_truth.centroids[i].x,
                   ground_truth.centroids[i].y,
                   ground_truth.centroids[i].flux);
        }
    }
    printf("\n");

    /* ========================================================= */
    /*  STEP 5: Build star pair database for identification       */
    /* ========================================================= */
    printf("[5] STAR PAIR DATABASE\n");
    float fov_rad = cam.fov_x * (float)(M_PI / 180.0) * 1.5f;  /* slightly wider */
    star_pair_db_build(&pair_db, &catalog, 5.0f, fov_rad);
    printf("    Pairs built: %u (mag <= 5.0, FOV %.1f deg)\n\n",
           pair_db.count, fov_rad * 57.2957795f);

    /* ========================================================= */
    /*  STEP 6: Lost-in-space star identification                 */
    /* ========================================================= */
    printf("[6] STAR IDENTIFICATION (Lost-in-Space)\n");
    StarIDResult id_result;
    int nmatched = star_id_lost_in_space(&ground_truth, &cam, &pair_db, &catalog, &id_result);
    printf("    Stars identified: %d\n", nmatched);
    for (uint16_t i = 0; i < id_result.count; i++) {
        uint16_t ci = id_result.matches[i].centroid_idx;
        uint16_t si = id_result.matches[i].catalog_idx;
        printf("    Centroid %u -> HIP %u (RA=%.3f, Dec=%+.3f, Mag=%.2f)\n",
               ci, catalog.stars[si].hip_id,
               catalog.stars[si].ra, catalog.stars[si].dec, catalog.stars[si].mag);
    }
    printf("\n");

    /* ========================================================= */
    /*  STEP 7: Attitude from ground truth (no noise)             */
    /* ========================================================= */
    printf("[7] ATTITUDE (Ground Truth - No Noise)\n");
    if (id_result.count >= 2) {
        MatchList ml_clean;
        match_list_init(&ml_clean);

        for (uint16_t i = 0; i < id_result.count; i++) {
            uint16_t ci = id_result.matches[i].centroid_idx;
            uint16_t si = id_result.matches[i].catalog_idx;

            float ux, uy, uz;
            pixel_to_unit_vector(&cam,
                                 ground_truth.centroids[ci].x,
                                 ground_truth.centroids[ci].y,
                                 &ux, &uy, &uz);

            match_list_add(&ml_clean,
                           ux, uy, uz,
                           catalog.stars[si].x_unit,
                           catalog.stars[si].y_unit,
                           catalog.stars[si].z_unit,
                           1.0f);
        }

        Quaternion q_clean;
        attitude_quest(&ml_clean, &q_clean);
        float err_clean = quaternion_error_angle(&true_q, &q_clean);

        DCM dcm_clean;
        quaternion_to_dcm(&q_clean, &dcm_clean);
        float r_c, p_c, y_c;
        dcm_to_euler(&dcm_clean, &r_c, &p_c, &y_c);

        printf("    QUEST Euler: Roll=%+.4f  Pitch=%+.4f  Yaw=%+.4f\n", r_c, p_c, y_c);
        printf("    Attitude error: %.6f deg (%.2f arcsec)\n\n", err_clean, err_clean * 3600.0f);
    }

    /* ========================================================= */
    /*  STEP 8: Noisy simulation - NOISE VALUES FROM USER         */
    /* ========================================================= */
    printf("[8] NOISE SIMULATION\n");
    printf("    Enter noise sigma values (pixels) to test.\n");
    printf("    Format: comma-separated, e.g. 0.1,0.3,0.5,1.0\n");
    printf("    > ");

    char input[256];
    if (fgets(input, sizeof(input), stdin) == NULL) {
        printf("    No input, using defaults: 0.1, 0.5, 1.0\n");
        strcpy(input, "0.1,0.5,1.0");
    }

    /* Parse noise values */
    float noise_vals[16];
    int n_noise = 0;
    char *tok = input;
    while (*tok && n_noise < 16) {
        float val = (float)atof(tok);
        if (val > 0.0f) {
            noise_vals[n_noise++] = val;
        }
        while (*tok && *tok != ',') tok++;
        if (*tok == ',') tok++;
    }

    if (n_noise == 0) {
        noise_vals[0] = 0.1f;
        noise_vals[1] = 0.5f;
        noise_vals[2] = 1.0f;
        n_noise = 3;
    }

    printf("\n    %-12s  %-14s  %-10s  %-14s  %-14s\n",
           "Noise(px)", "Centroid RMS", "Stars ID", "Att.Err(deg)", "Att.Err(arcsec)");
    printf("    %-12s  %-14s  %-10s  %-14s  %-14s\n",
           "----------", "-----------", "--------", "-----------", "--------------");

    for (int ni = 0; ni < n_noise; ni++) {
        float sigma = noise_vals[ni];

        /* Copy ground truth and add noise */
        CentroidList noisy;
        memcpy(&noisy, &ground_truth, sizeof(CentroidList));
        add_centroid_noise(&noisy, sigma, 42 + (uint32_t)(ni * 1000));

        /* Compute centroid RMS error */
        float rms = compute_centroid_rms_error(&noisy, &ground_truth);

        /* Star ID on noisy data */
        StarIDResult noisy_id;
        int n_id = star_id_lost_in_space(&noisy, &cam, &pair_db, &catalog, &noisy_id);

        /* Attitude from noisy centroids */
        float att_err_deg = -1.0f;
        float att_err_arcsec = -1.0f;
        SimulationResult sim;
        memset(&sim, 0, sizeof(sim));
        sim.noise_sigma = sigma;
        sim.centroid_rms = rms;
        sim.stars_identified = n_id;

        if (noisy_id.count >= 2) {
            MatchList ml_noisy;
            match_list_init(&ml_noisy);

            for (uint16_t i = 0; i < noisy_id.count; i++) {
                uint16_t ci = noisy_id.matches[i].centroid_idx;
                uint16_t si = noisy_id.matches[i].catalog_idx;

                float ux, uy, uz;
                pixel_to_unit_vector(&cam,
                                     noisy.centroids[ci].x,
                                     noisy.centroids[ci].y,
                                     &ux, &uy, &uz);

                match_list_add(&ml_noisy,
                               ux, uy, uz,
                               catalog.stars[si].x_unit,
                               catalog.stars[si].y_unit,
                               catalog.stars[si].z_unit,
                               1.0f);
            }

            Quaternion q_noisy;
            attitude_quest(&ml_noisy, &q_noisy);
            att_err_deg = quaternion_error_angle(&true_q, &q_noisy);
            att_err_arcsec = att_err_deg * 3600.0f;

            DCM dcm_noisy;
            quaternion_to_dcm(&q_noisy, &dcm_noisy);
            float rn, pn, yn;
            dcm_to_euler(&dcm_noisy, &rn, &pn, &yn);
            sim.roll_error = rn - true_roll;
            sim.pitch_error = pn - true_pitch;
            sim.yaw_error = yn - true_yaw;
        }

        sim.attitude_error_deg = att_err_deg;
        sim.attitude_error_arcsec = att_err_arcsec;

        if (att_err_deg >= 0.0f) {
            printf("    %-12.4f  %-14.4f  %-10d  %-14.6f  %-14.2f\n",
                   sigma, rms, n_id, att_err_deg, att_err_arcsec);
        } else {
            printf("    %-12.4f  %-14.4f  %-10d  FAILED          FAILED\n",
                   sigma, rms, n_id);
        }
    }

    /* ========================================================= */
    /*  STEP 9: Detailed error report for last noise value        */
    /* ========================================================= */
    if (n_noise > 0) {
        float sigma = noise_vals[n_noise - 1];
        CentroidList noisy_last;
        memcpy(&noisy_last, &ground_truth, sizeof(CentroidList));
        add_centroid_noise(&noisy_last, sigma, 42 + (uint32_t)((n_noise - 1) * 1000));

        printf("\n[9] DETAILED ERROR ANALYSIS (sigma = %.4f px)\n", sigma);

        float errors_x[MAX_CENTROIDS_PER_FRAME];
        float errors_y[MAX_CENTROIDS_PER_FRAME];
        uint16_t nerr;
        compute_centroid_errors(&noisy_last, &ground_truth, errors_x, errors_y, &nerr);

        printf("    %-4s  %-12s  %-12s  %-12s  %-12s  %-12s\n",
               "ID", "GT_X", "GT_Y", "Noisy_X", "Noisy_Y", "Error(px)");
        int show = nerr > 8 ? 8 : nerr;
        for (int i = 0; i < show; i++) {
            float err = sqrtf(errors_x[i]*errors_x[i] + errors_y[i]*errors_y[i]);
            printf("    %-4u  %-12.4f  %-12.4f  %-12.4f  %-12.4f  %-12.4f\n",
                   i,
                   ground_truth.centroids[i].x, ground_truth.centroids[i].y,
                   noisy_last.centroids[i].x, noisy_last.centroids[i].y,
                   err);
        }
    }

    printf("\n*** Pipeline Complete ***\n");
    return 0;
}
