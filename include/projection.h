#ifndef PROJECTION_H
#define PROJECTION_H

#include "star_tracker.h"
#include "catalog.h"
#include "centroid.h"
#include "attitude.h"

/*
 * Project catalog stars onto the image plane given a camera and attitude.
 * Only stars within the FOV are projected.
 *
 * cam       - camera parameters
 * dcm       - attitude DCM (inertial to body rotation)
 * catalog   - star catalog
 * mag_limit - only project stars brighter than this
 * out       - output centroid list (ground truth, sub-pixel accuracy)
 *
 * Returns number of stars projected.
 */
int project_stars_to_image(const CameraParams *cam,
                           const DCM *dcm,
                           const StarCatalog *catalog,
                           float mag_limit,
                           CentroidList *out);

/*
 * Convert a single inertial unit vector to pixel coordinates
 * given camera and attitude.
 * Returns 0 if the star is within the image, -1 if outside.
 */
int unit_vector_to_pixel(const CameraParams *cam,
                         const DCM *dcm,
                         float ref_x, float ref_y, float ref_z,
                         float *px, float *py);

/*
 * Flux model: approximate flux from magnitude
 * Uses: flux = base_flux * 10^(-0.4 * mag)
 */
float magnitude_to_flux(float mag, float base_flux);

#endif /* PROJECTION_H */
