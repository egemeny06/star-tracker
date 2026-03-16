#include <stdio.h>
#include <math.h>
#include "centroid.h"
#include "catalog.h"
#include "star_tracker.h"
#include "attitude.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/* Sample Hipparcos stars for testing (brightest stars) */
static void load_sample_catalog(StarCatalog *catalog)
{
    catalog_init(catalog);

    /* HIP ID,  RA(deg),   Dec(deg),  Mag  */
    catalog_add_star(catalog, 32349,  101.287, -16.716, -1.46f); /* Sirius */
    catalog_add_star(catalog, 30438,   95.988, -52.696, -0.74f); /* Canopus */
    catalog_add_star(catalog, 69673,  213.915,  19.182, -0.05f); /* Arcturus */
    catalog_add_star(catalog, 71683,  219.902, -60.834, -0.01f); /* Alpha Centauri A */
    catalog_add_star(catalog, 91262,  279.235,  38.784,  0.03f); /* Vega */
    catalog_add_star(catalog, 24436,   78.634,  -8.202,  0.12f); /* Rigel */
    catalog_add_star(catalog, 37279,  114.827,   5.225,  0.38f); /* Procyon */
    catalog_add_star(catalog, 24608,   79.172,  45.998,  0.08f); /* Capella */
    catalog_add_star(catalog, 27989,   88.793,   7.407,  0.50f); /* Betelgeuse */
    catalog_add_star(catalog, 7588,   24.429, -57.237,  0.46f); /* Achernar */
}

/* Build a known rotation matrix from Euler angles (deg) for testing */
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

/* Rotate a reference vector by DCM to simulate a body-frame observation */
static void rotate_vec(const DCM *dcm, const float *vin, float *vout)
{
    for (int i = 0; i < 3; i++) {
        vout[i] = dcm->m[i][0] * vin[0]
                + dcm->m[i][1] * vin[1]
                + dcm->m[i][2] * vin[2];
    }
}

int main(void)
{
    printf("*** Star Tracker - Centroid & Attitude Demo ***\n\n");

    /* --- 1. Camera Setup --- */
    CameraParams cam;
    camera_init(&cam, 50.0f, 0.0055f, 1024, 1024);
    printf("Camera: %.1f mm focal, %.4f mm/px, %ux%u\n",
           cam.focal_length, cam.pixel_size, cam.width, cam.height);
    printf("FOV: %.2f x %.2f degrees\n\n", cam.fov_x, cam.fov_y);

    /* --- 2. Load Hipparcos Sample Catalog --- */
    StarCatalog catalog;
    load_sample_catalog(&catalog);
    printf("Catalog loaded: %u stars\n\n", catalog.count);

    /* --- 3. Centroid extraction demo --- */
    CentroidList frame;
    centroid_list_init(&frame, 1);
    centroid_list_add(&frame, 512.3f, 510.7f, 45000.0f, 0.12f, 0.11f);
    centroid_list_add(&frame, 234.1f, 678.9f, 32000.0f, 0.15f, 0.14f);
    centroid_list_add(&frame, 789.5f, 123.4f, 28500.0f, 0.18f, 0.17f);
    centroid_list_add(&frame, 100.2f, 900.1f, 15000.0f, 0.25f, 0.22f);
    centroid_list_add(&frame, 950.8f, 450.3f, 12000.0f, 0.30f, 0.28f);
    centroid_list_add(&frame, 400.0f, 300.0f,  8500.0f, 0.35f, 0.33f);

    centroid_list_sort_by_flux(&frame);
    printf("--- Centroids (sorted by flux) ---\n");
    centroid_list_print(&frame);

    /* ============================================================== */
    /*  ATTITUDE DETERMINATION TEST                                   */
    /*  We define a known attitude, rotate catalog stars to simulate   */
    /*  body-frame observations, then recover the attitude.           */
    /* ============================================================== */

    printf("========================================\n");
    printf("  ATTITUDE DETERMINATION TEST\n");
    printf("========================================\n\n");

    /* --- 4. Define TRUE attitude (what we want to recover) --- */
    float true_roll  = 15.0f;
    float true_pitch = -25.0f;
    float true_yaw   = 40.0f;

    DCM true_dcm;
    euler_to_dcm(true_roll, true_pitch, true_yaw, &true_dcm);

    printf("TRUE attitude (Euler angles):\n");
    printf("  Roll:  %+.2f deg\n", true_roll);
    printf("  Pitch: %+.2f deg\n", true_pitch);
    printf("  Yaw:   %+.2f deg\n\n", true_yaw);

    printf("TRUE DCM:\n");
    dcm_print(&true_dcm);

    Quaternion true_q;
    dcm_to_quaternion(&true_dcm, &true_q);
    printf("TRUE Quaternion: ");
    quaternion_print(&true_q);
    printf("\n");

    /* --- 5. Simulate star observations --- */
    /* Pick 5 catalog stars, get their inertial unit vectors,
       rotate them by true_dcm to get "body-frame" observations */
    printf("--- Simulating %d star observations ---\n", 5);

    MatchList matches;
    match_list_init(&matches);

    for (int i = 0; i < 5 && i < catalog.count; i++) {
        float ref[3] = {
            catalog.stars[i].x_unit,
            catalog.stars[i].y_unit,
            catalog.stars[i].z_unit
        };
        float body[3];
        rotate_vec(&true_dcm, ref, body);

        printf("  Star HIP %5u: ref=(%+.4f, %+.4f, %+.4f) -> body=(%+.4f, %+.4f, %+.4f)\n",
               catalog.stars[i].hip_id,
               ref[0], ref[1], ref[2],
               body[0], body[1], body[2]);

        match_list_add(&matches,
                       body[0], body[1], body[2],
                       ref[0], ref[1], ref[2],
                       1.0f);
    }
    printf("\n");

    /* --- 6. TRIAD Algorithm --- */
    printf("--- TRIAD Algorithm (2 vectors) ---\n");
    DCM triad_dcm;
    if (attitude_triad(&matches, &triad_dcm) == 0) {
        dcm_print(&triad_dcm);

        Quaternion triad_q;
        dcm_to_quaternion(&triad_dcm, &triad_q);
        printf("TRIAD Quaternion: ");
        quaternion_print(&triad_q);

        float roll, pitch, yaw;
        dcm_to_euler(&triad_dcm, &roll, &pitch, &yaw);
        printf("TRIAD Euler: Roll=%+.4f  Pitch=%+.4f  Yaw=%+.4f\n", roll, pitch, yaw);

        float err = quaternion_error_angle(&true_q, &triad_q);
        printf("TRIAD error vs truth: %.6f deg\n", err);
    } else {
        printf("TRIAD failed!\n");
    }
    printf("\n");

    /* --- 7. QUEST Algorithm --- */
    printf("--- QUEST Algorithm (%u vectors) ---\n", matches.count);
    Quaternion quest_q;
    if (attitude_quest(&matches, &quest_q) == 0) {
        printf("QUEST Quaternion: ");
        quaternion_print(&quest_q);

        DCM quest_dcm;
        quaternion_to_dcm(&quest_q, &quest_dcm);
        dcm_print(&quest_dcm);

        float roll, pitch, yaw;
        dcm_to_euler(&quest_dcm, &roll, &pitch, &yaw);
        printf("QUEST Euler: Roll=%+.4f  Pitch=%+.4f  Yaw=%+.4f\n", roll, pitch, yaw);

        float err = quaternion_error_angle(&true_q, &quest_q);
        printf("QUEST error vs truth: %.6f deg\n", err);
    } else {
        printf("QUEST failed!\n");
    }
    printf("\n");

    /* --- 8. Comparison summary --- */
    printf("========================================\n");
    printf("  SUMMARY\n");
    printf("========================================\n");
    printf("True attitude:  Roll=%+.2f  Pitch=%+.2f  Yaw=%+.2f\n",
           true_roll, true_pitch, true_yaw);

    {
        float r, p, y;
        dcm_to_euler(&triad_dcm, &r, &p, &y);
        printf("TRIAD result:   Roll=%+.4f  Pitch=%+.4f  Yaw=%+.4f\n", r, p, y);
    }
    {
        DCM qd;
        quaternion_to_dcm(&quest_q, &qd);
        float r, p, y;
        dcm_to_euler(&qd, &r, &p, &y);
        printf("QUEST result:   Roll=%+.4f  Pitch=%+.4f  Yaw=%+.4f\n", r, p, y);
    }

    printf("\n*** Done ***\n");
    return 0;
}
