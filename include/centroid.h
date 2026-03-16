#ifndef CENTROID_H
#define CENTROID_H

#include <stdint.h>
#include <stddef.h>

#define MAX_CENTROIDS_PER_FRAME 256

typedef struct {
    float x;          /* centroid x position (pixels) */
    float y;          /* centroid y position (pixels) */
    float flux;       /* integrated brightness */
    float x_err;      /* centroid x error (pixels) */
    float y_err;      /* centroid y error (pixels) */
    uint16_t id;      /* centroid index within frame */
} Centroid;

typedef struct {
    Centroid centroids[MAX_CENTROIDS_PER_FRAME];
    uint16_t count;
    uint32_t frame_id;
} CentroidList;

/* Initialize a centroid list for a given frame */
void centroid_list_init(CentroidList *list, uint32_t frame_id);

/* Add a centroid to the list. Returns 0 on success, -1 if full. */
int centroid_list_add(CentroidList *list, float x, float y, float flux,
                      float x_err, float y_err);

/* Sort centroids by flux (brightest first) */
void centroid_list_sort_by_flux(CentroidList *list);

/* Print centroid list to stdout (debug) */
void centroid_list_print(const CentroidList *list);

#endif /* CENTROID_H */
