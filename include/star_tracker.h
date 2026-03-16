#ifndef STAR_TRACKER_H
#define STAR_TRACKER_H

#include "centroid.h"
#include "catalog.h"

/* Camera / sensor parameters */
typedef struct {
    float focal_length;    /* mm */
    float pixel_size;      /* mm/pixel */
    uint16_t width;        /* sensor width (pixels) */
    uint16_t height;       /* sensor height (pixels) */
    float fov_x;           /* field of view x (degrees) */
    float fov_y;           /* field of view y (degrees) */
} CameraParams;

/* Initialize camera parameters and compute FOV */
void camera_init(CameraParams *cam, float focal_length_mm,
                 float pixel_size_mm, uint16_t width, uint16_t height);

/* Convert pixel coordinates to unit vector in camera frame */
void pixel_to_unit_vector(const CameraParams *cam, float px, float py,
                          float *ux, float *uy, float *uz);

/* Compute angular distance between two unit vectors (radians) */
float angular_distance(float ux1, float uy1, float uz1,
                       float ux2, float uy2, float uz2);

#endif /* STAR_TRACKER_H */
