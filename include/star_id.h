#ifndef STAR_ID_H
#define STAR_ID_H

#include "centroid.h"
#include "catalog.h"
#include "star_tracker.h"
#include "attitude.h"

#define MAX_STAR_PAIRS 16384

/* A pair of catalog stars with their precomputed angular distance */
typedef struct {
    uint16_t idx_a;        /* catalog index of star A */
    uint16_t idx_b;        /* catalog index of star B */
    float ang_dist;        /* angular distance (radians) */
} StarPair;

/* Star pair database sorted by angular distance for binary search */
typedef struct {
    StarPair pairs[MAX_STAR_PAIRS];
    uint16_t count;
} StarPairDB;

/* A candidate match: centroid index → catalog index */
typedef struct {
    uint16_t centroid_idx;
    uint16_t catalog_idx;
    float score;           /* match confidence */
} StarIDMatch;

#define MAX_ID_MATCHES 64

typedef struct {
    StarIDMatch matches[MAX_ID_MATCHES];
    uint16_t count;
} StarIDResult;

/*
 * Build a star pair database from the catalog.
 * Only includes stars brighter than mag_limit and pairs within max_fov (rad).
 */
void star_pair_db_build(StarPairDB *db,
                        const StarCatalog *catalog,
                        float mag_limit,
                        float max_fov_rad);

/*
 * Binary search in star pair DB for pairs matching a given angular distance
 * within tolerance. Returns number of matches found.
 * out_indices: indices into db->pairs array, max_out: max results
 */
int star_pair_db_search(const StarPairDB *db,
                        float target_dist,
                        float tolerance,
                        uint16_t *out_indices,
                        int max_out);

/*
 * Lost-in-space star identification.
 * Takes observed centroids, camera params, and star pair DB.
 * Returns matched centroid-to-catalog pairs.
 */
int star_id_lost_in_space(const CentroidList *centroids,
                          const CameraParams *cam,
                          const StarPairDB *db,
                          const StarCatalog *catalog,
                          StarIDResult *result);

#endif /* STAR_ID_H */
