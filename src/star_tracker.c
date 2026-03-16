#include "star_tracker.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define RAD2DEG(r) ((r) * (float)(180.0 / M_PI))

void camera_init(CameraParams *cam, float focal_length_mm,
                 float pixel_size_mm, uint16_t width, uint16_t height)
{
    cam->focal_length = focal_length_mm;
    cam->pixel_size = pixel_size_mm;
    cam->width = width;
    cam->height = height;

    /* Compute field of view */
    float sensor_w = width * pixel_size_mm;
    float sensor_h = height * pixel_size_mm;
    cam->fov_x = 2.0f * RAD2DEG(atanf(sensor_w / (2.0f * focal_length_mm)));
    cam->fov_y = 2.0f * RAD2DEG(atanf(sensor_h / (2.0f * focal_length_mm)));
}

void pixel_to_unit_vector(const CameraParams *cam, float px, float py,
                          float *ux, float *uy, float *uz)
{
    /* Convert pixel to camera-frame coordinates (origin at sensor center) */
    float cx = (px - cam->width / 2.0f) * cam->pixel_size;
    float cy = (py - cam->height / 2.0f) * cam->pixel_size;
    float f = cam->focal_length;

    /* Unit vector: boresight along Z */
    float norm = sqrtf(cx * cx + cy * cy + f * f);
    *ux = cx / norm;
    *uy = cy / norm;
    *uz = f / norm;
}

float angular_distance(float ux1, float uy1, float uz1,
                       float ux2, float uy2, float uz2)
{
    float dot = ux1 * ux2 + uy1 * uy2 + uz1 * uz2;

    /* Clamp for numerical safety */
    if (dot > 1.0f) dot = 1.0f;
    if (dot < -1.0f) dot = -1.0f;

    return acosf(dot);
}
