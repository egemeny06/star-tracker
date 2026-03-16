#include "catalog.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define DEG2RAD(d) ((d) * (float)(M_PI / 180.0))

void catalog_init(StarCatalog *catalog)
{
    memset(catalog, 0, sizeof(StarCatalog));
    catalog->count = 0;
}

int catalog_add_star(StarCatalog *catalog, uint32_t hip_id,
                     float ra_deg, float dec_deg, float mag)
{
    if (catalog->count >= CATALOG_MAX_STARS)
        return -1;

    CatalogStar *s = &catalog->stars[catalog->count];
    s->hip_id = hip_id;
    s->ra = ra_deg;
    s->dec = dec_deg;
    s->mag = mag;

    /* Precompute unit vector */
    float ra_rad = DEG2RAD(ra_deg);
    float dec_rad = DEG2RAD(dec_deg);
    s->x_unit = cosf(dec_rad) * cosf(ra_rad);
    s->y_unit = cosf(dec_rad) * sinf(ra_rad);
    s->z_unit = sinf(dec_rad);

    catalog->count++;
    return 0;
}

int catalog_load_csv(StarCatalog *catalog, const char *filepath)
{
    FILE *fp = fopen(filepath, "r");
    if (!fp)
        return -1;

    char line[256];
    int loaded = 0;

    /* Skip header line */
    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return -1;
    }

    while (fgets(line, sizeof(line), fp) != NULL) {
        uint32_t hip_id;
        float ra, dec, mag;

        if (sscanf(line, "%u,%f,%f,%f", &hip_id, &ra, &dec, &mag) == 4) {
            if (catalog_add_star(catalog, hip_id, ra, dec, mag) == 0)
                loaded++;
        }
    }

    fclose(fp);
    return loaded;
}

uint16_t catalog_filter_by_magnitude(const StarCatalog *catalog,
                                     float mag_limit,
                                     CatalogStar *out, uint16_t max_out)
{
    uint16_t count = 0;

    for (uint16_t i = 0; i < catalog->count && count < max_out; i++) {
        if (catalog->stars[i].mag <= mag_limit) {
            out[count] = catalog->stars[i];
            count++;
        }
    }
    return count;
}

void catalog_compute_unit_vectors(StarCatalog *catalog)
{
    for (uint16_t i = 0; i < catalog->count; i++) {
        CatalogStar *s = &catalog->stars[i];
        float ra_rad = DEG2RAD(s->ra);
        float dec_rad = DEG2RAD(s->dec);
        s->x_unit = cosf(dec_rad) * cosf(ra_rad);
        s->y_unit = cosf(dec_rad) * sinf(ra_rad);
        s->z_unit = sinf(dec_rad);
    }
}
