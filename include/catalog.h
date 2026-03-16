#ifndef CATALOG_H
#define CATALOG_H

#include <stdint.h>
#include <stddef.h>

#define CATALOG_MAX_STARS 9110  /* Hipparcos bright star count */

typedef struct {
    uint32_t hip_id;       /* Hipparcos catalog number */
    float ra;              /* Right ascension (degrees) */
    float dec;             /* Declination (degrees) */
    float mag;             /* Visual magnitude */
    float x_unit;          /* Unit vector x (precomputed) */
    float y_unit;          /* Unit vector y (precomputed) */
    float z_unit;          /* Unit vector z (precomputed) */
} CatalogStar;

typedef struct {
    CatalogStar stars[CATALOG_MAX_STARS];
    uint16_t count;
} StarCatalog;

/* Initialize catalog (zeroes out) */
void catalog_init(StarCatalog *catalog);

/* Load catalog from CSV file. Returns number of stars loaded, -1 on error. */
int catalog_load_csv(StarCatalog *catalog, const char *filepath);

/* Add a star manually (for testing or embedded ROM catalog) */
int catalog_add_star(StarCatalog *catalog, uint32_t hip_id,
                     float ra_deg, float dec_deg, float mag);

/* Get stars brighter than a magnitude limit */
uint16_t catalog_filter_by_magnitude(const StarCatalog *catalog,
                                     float mag_limit,
                                     CatalogStar *out, uint16_t max_out);

/* Precompute unit vectors from RA/Dec for all stars */
void catalog_compute_unit_vectors(StarCatalog *catalog);

#endif /* CATALOG_H */
