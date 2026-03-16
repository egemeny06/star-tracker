#include "centroid.h"
#include <stdio.h>
#include <string.h>

void centroid_list_init(CentroidList *list, uint32_t frame_id)
{
    memset(list, 0, sizeof(CentroidList));
    list->frame_id = frame_id;
    list->count = 0;
}

int centroid_list_add(CentroidList *list, float x, float y, float flux,
                      float x_err, float y_err)
{
    if (list->count >= MAX_CENTROIDS_PER_FRAME)
        return -1;

    Centroid *c = &list->centroids[list->count];
    c->x = x;
    c->y = y;
    c->flux = flux;
    c->x_err = x_err;
    c->y_err = y_err;
    c->id = list->count;
    list->count++;
    return 0;
}

void centroid_list_sort_by_flux(CentroidList *list)
{
    /* Simple insertion sort - embedded friendly, no malloc */
    uint16_t i, j;
    Centroid temp;

    for (i = 1; i < list->count; i++) {
        temp = list->centroids[i];
        j = i;
        while (j > 0 && list->centroids[j - 1].flux < temp.flux) {
            list->centroids[j] = list->centroids[j - 1];
            j--;
        }
        list->centroids[j] = temp;
    }
}

void centroid_list_print(const CentroidList *list)
{
    printf("=== Frame %u | Centroid Count: %u ===\n",
           list->frame_id, list->count);
    printf("%-4s  %-10s  %-10s  %-12s  %-10s  %-10s\n",
           "ID", "X (px)", "Y (px)", "Flux", "X_err", "Y_err");
    printf("----  ----------  ----------  ------------  ----------  ----------\n");

    for (uint16_t i = 0; i < list->count; i++) {
        const Centroid *c = &list->centroids[i];
        printf("%-4u  %-10.3f  %-10.3f  %-12.2f  %-10.4f  %-10.4f\n",
               c->id, c->x, c->y, c->flux, c->x_err, c->y_err);
    }
    printf("\n");
}
