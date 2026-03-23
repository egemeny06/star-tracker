#include "star_id.h"
#include "star_tracker.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* ------------------------------------------------------------------ */
/*  Insertion sort for star pairs by angular distance                  */
/* ------------------------------------------------------------------ */
static void sort_pairs(StarPair *pairs, uint16_t count)
{
    for (uint16_t i = 1; i < count; i++) {
        StarPair temp = pairs[i];
        uint16_t j = i;
        while (j > 0 && pairs[j - 1].ang_dist > temp.ang_dist) {
            pairs[j] = pairs[j - 1];
            j--;
        }
        pairs[j] = temp;
    }
}

/* ------------------------------------------------------------------ */
/*  Build star pair database                                           */
/* ------------------------------------------------------------------ */
void star_pair_db_build(StarPairDB *db,
                        const StarCatalog *catalog,
                        float mag_limit,
                        float max_fov_rad)
{
    memset(db, 0, sizeof(StarPairDB));
    db->count = 0;

    /* First, collect indices of bright stars */
    uint16_t bright_idx[1024];
    uint16_t nbright = 0;

    for (uint16_t i = 0; i < catalog->count && nbright < 1024; i++) {
        if (catalog->stars[i].mag <= mag_limit) {
            bright_idx[nbright++] = i;
        }
    }

    /* Build pairs within FOV limit */
    for (uint16_t i = 0; i < nbright && db->count < MAX_STAR_PAIRS; i++) {
        const CatalogStar *a = &catalog->stars[bright_idx[i]];
        for (uint16_t j = i + 1; j < nbright && db->count < MAX_STAR_PAIRS; j++) {
            const CatalogStar *b = &catalog->stars[bright_idx[j]];

            float dist = angular_distance(a->x_unit, a->y_unit, a->z_unit,
                                          b->x_unit, b->y_unit, b->z_unit);

            if (dist <= max_fov_rad && dist > 0.001f) {
                StarPair *p = &db->pairs[db->count];
                p->idx_a = bright_idx[i];
                p->idx_b = bright_idx[j];
                p->ang_dist = dist;
                db->count++;
            }
        }
    }

    /* Sort by angular distance for binary search */
    sort_pairs(db->pairs, db->count);
}

/* ------------------------------------------------------------------ */
/*  Binary search for angular distance match                           */
/* ------------------------------------------------------------------ */
int star_pair_db_search(const StarPairDB *db,
                        float target_dist,
                        float tolerance,
                        uint16_t *out_indices,
                        int max_out)
{
    float lo = target_dist - tolerance;
    float hi = target_dist + tolerance;
    int found = 0;

    /* Binary search for lower bound */
    int left = 0, right = db->count - 1;
    int start = db->count;

    while (left <= right) {
        int mid = (left + right) / 2;
        if (db->pairs[mid].ang_dist >= lo) {
            start = mid;
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }

    /* Collect all matches in range [lo, hi] */
    for (int i = start; i < db->count && found < max_out; i++) {
        if (db->pairs[i].ang_dist > hi)
            break;
        out_indices[found++] = (uint16_t)i;
    }

    return found;
}

/* ------------------------------------------------------------------ */
/*  Lost-in-space identification                                       */
/*  Strategy: take brightest centroid pair, compute angular distance,   */
/*  search pair DB, verify with third star.                            */
/* ------------------------------------------------------------------ */
int star_id_lost_in_space(const CentroidList *centroids,
                          const CameraParams *cam,
                          const StarPairDB *db,
                          const StarCatalog *catalog,
                          StarIDResult *result)
{
    memset(result, 0, sizeof(StarIDResult));

    if (centroids->count < 3)
        return -1;

    /* Convert top centroids to unit vectors */
    int nuse = centroids->count;
    if (nuse > 10) nuse = 10;  /* limit for speed */

    float uvecs[10][3];
    for (int i = 0; i < nuse; i++) {
        pixel_to_unit_vector(cam,
                             centroids->centroids[i].x,
                             centroids->centroids[i].y,
                             &uvecs[i][0], &uvecs[i][1], &uvecs[i][2]);
    }

    /* Compute observed angular distances between pairs */
    float obs_dist_01 = angular_distance(uvecs[0][0], uvecs[0][1], uvecs[0][2],
                                         uvecs[1][0], uvecs[1][1], uvecs[1][2]);
    float obs_dist_02 = angular_distance(uvecs[0][0], uvecs[0][1], uvecs[0][2],
                                         uvecs[2][0], uvecs[2][1], uvecs[2][2]);
    float obs_dist_12 = angular_distance(uvecs[1][0], uvecs[1][1], uvecs[1][2],
                                         uvecs[2][0], uvecs[2][1], uvecs[2][2]);

    /* Search pair DB for distance 0-1 */
    float tol = 0.0005f;  /* ~0.03 degrees tolerance */
    uint16_t match_01[64];
    int n01 = star_pair_db_search(db, obs_dist_01, tol, match_01, 64);

    if (n01 == 0) {
        tol = 0.002f;  /* widen tolerance */
        n01 = star_pair_db_search(db, obs_dist_01, tol, match_01, 64);
    }

    /* For each candidate pair matching 0-1, verify with star 2 */
    for (int k = 0; k < n01; k++) {
        const StarPair *p01 = &db->pairs[match_01[k]];
        uint16_t cat_a = p01->idx_a;
        uint16_t cat_b = p01->idx_b;

        const CatalogStar *sa = &catalog->stars[cat_a];
        const CatalogStar *sb = &catalog->stars[cat_b];

        /* Check: does star 2 match angular distances to both a and b? */
        /* Try both orientations: centroid0=a,centroid1=b and centroid0=b,centroid1=a */
        for (int flip = 0; flip < 2; flip++) {
            uint16_t c0 = (flip == 0) ? cat_a : cat_b;
            uint16_t c1 = (flip == 0) ? cat_b : cat_a;
            const CatalogStar *s0 = &catalog->stars[c0];
            const CatalogStar *s1 = &catalog->stars[c1];

            /* Search for star 2: find catalog star matching dist to both s0 and s1 */
            uint16_t match_02[32];
            int n02 = star_pair_db_search(db, obs_dist_02, tol, match_02, 32);

            for (int m = 0; m < n02; m++) {
                const StarPair *p02 = &db->pairs[match_02[m]];

                /* Which end matches s0? */
                uint16_t cat_2 = 0xFFFF;
                if (p02->idx_a == c0)
                    cat_2 = p02->idx_b;
                else if (p02->idx_b == c0)
                    cat_2 = p02->idx_a;
                else
                    continue;

                /* Verify distance 1-2 */
                const CatalogStar *s2 = &catalog->stars[cat_2];
                float cat_dist_12 = angular_distance(
                    s1->x_unit, s1->y_unit, s1->z_unit,
                    s2->x_unit, s2->y_unit, s2->z_unit);

                if (fabsf(cat_dist_12 - obs_dist_12) < tol) {
                    /* Found a valid triple match! */
                    (void)sa; (void)sb; (void)s0;

                    if (result->count < MAX_ID_MATCHES) {
                        result->matches[result->count].centroid_idx = 0;
                        result->matches[result->count].catalog_idx = c0;
                        result->matches[result->count].score = 1.0f;
                        result->count++;
                    }
                    if (result->count < MAX_ID_MATCHES) {
                        result->matches[result->count].centroid_idx = 1;
                        result->matches[result->count].catalog_idx = c1;
                        result->matches[result->count].score = 1.0f;
                        result->count++;
                    }
                    if (result->count < MAX_ID_MATCHES) {
                        result->matches[result->count].centroid_idx = 2;
                        result->matches[result->count].catalog_idx = cat_2;
                        result->matches[result->count].score = 1.0f;
                        result->count++;
                    }

                    /* Try to match remaining centroids */
                    for (int ci = 3; ci < nuse && result->count < MAX_ID_MATCHES; ci++) {
                        float best_err = 1e6f;
                        uint16_t best_cat = 0xFFFF;

                        for (uint16_t si = 0; si < catalog->count; si++) {
                            if (catalog->stars[si].mag > 6.0f) continue;

                            const CatalogStar *sc = &catalog->stars[si];
                            float d0 = angular_distance(
                                s0->x_unit, s0->y_unit, s0->z_unit,
                                sc->x_unit, sc->y_unit, sc->z_unit);
                            float obs_d0 = angular_distance(
                                uvecs[0][0], uvecs[0][1], uvecs[0][2],
                                uvecs[ci][0], uvecs[ci][1], uvecs[ci][2]);

                            float err = fabsf(d0 - obs_d0);
                            if (err < tol && err < best_err) {
                                /* Double check with second reference */
                                float d1 = angular_distance(
                                    s1->x_unit, s1->y_unit, s1->z_unit,
                                    sc->x_unit, sc->y_unit, sc->z_unit);
                                float obs_d1 = angular_distance(
                                    uvecs[1][0], uvecs[1][1], uvecs[1][2],
                                    uvecs[ci][0], uvecs[ci][1], uvecs[ci][2]);

                                if (fabsf(d1 - obs_d1) < tol) {
                                    best_err = err;
                                    best_cat = si;
                                }
                            }
                        }

                        if (best_cat != 0xFFFF) {
                            result->matches[result->count].centroid_idx = (uint16_t)ci;
                            result->matches[result->count].catalog_idx = best_cat;
                            result->matches[result->count].score = 1.0f - best_err;
                            result->count++;
                        }
                    }

                    return (int)result->count;
                }
            }
        }
    }

    return 0; /* no match found */
}
