#include "projection.h"
#include <math.h>

float magnitude_to_flux(float mag, float base_flux)
{
    return base_flux * powf(10.0f, -0.4f * mag);
}

int unit_vector_to_pixel(const CameraParams *cam,
                         const DCM *dcm,
                         float ref_x, float ref_y, float ref_z,
                         float *px, float *py)
{
    /* Rotate inertial vector to body frame: b = DCM * r */
    float bx = dcm->m[0][0] * ref_x + dcm->m[0][1] * ref_y + dcm->m[0][2] * ref_z;
    float by = dcm->m[1][0] * ref_x + dcm->m[1][1] * ref_y + dcm->m[1][2] * ref_z;
    float bz = dcm->m[2][0] * ref_x + dcm->m[2][1] * ref_y + dcm->m[2][2] * ref_z;

    /* Star must be in front of camera (boresight = +Z) */
    if (bz <= 0.0f)
        return -1;

    /* Pinhole projection: body frame → image plane */
    float x_mm = cam->focal_length * bx / bz;
    float y_mm = cam->focal_length * by / bz;

    /* Convert mm to pixel coordinates (origin at top-left) */
    *px = x_mm / cam->pixel_size + cam->width / 2.0f;
    *py = y_mm / cam->pixel_size + cam->height / 2.0f;

    /* Check if within image bounds */
    if (*px < 0.0f || *px >= (float)cam->width ||
        *py < 0.0f || *py >= (float)cam->height)
        return -1;

    return 0;
}

int project_stars_to_image(const CameraParams *cam,
                           const DCM *dcm,
                           const StarCatalog *catalog,
                           float mag_limit,
                           CentroidList *out)
{
    centroid_list_init(out, 0);
    int count = 0;

    for (uint16_t i = 0; i < catalog->count; i++) {
        const CatalogStar *s = &catalog->stars[i];

        if (s->mag > mag_limit)
            continue;

        float px, py;
        if (unit_vector_to_pixel(cam, dcm,
                                 s->x_unit, s->y_unit, s->z_unit,
                                 &px, &py) == 0) {
            float flux = magnitude_to_flux(s->mag, 100000.0f);
            /* Ground truth: sub-pixel precision, zero error */
            if (centroid_list_add(out, px, py, flux, 0.0f, 0.0f) == 0) {
                count++;
            }
        }
    }

    return count;
}
