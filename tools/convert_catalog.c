/*
 * Convert HYG database CSV to star tracker catalog CSV.
 * Extracts only Hipparcos stars (hip > 0) and outputs:
 *   hip_id, ra_deg, dec_deg, mag
 *
 * Usage: convert_catalog <input_hyg.csv> <output.csv> [mag_limit]
 *   mag_limit defaults to 6.5 (naked eye stars)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 2048

/* Parse a quoted or unquoted CSV field, advance pointer */
static const char *next_field(const char *p, char *out, int maxlen)
{
    int i = 0;
    if (*p == '"') {
        p++;
        while (*p && !(*p == '"' && (*(p+1) == ',' || *(p+1) == '\0' || *(p+1) == '\n' || *(p+1) == '\r'))) {
            if (i < maxlen - 1) out[i++] = *p;
            p++;
        }
        if (*p == '"') p++;
        if (*p == ',') p++;
    } else {
        while (*p && *p != ',' && *p != '\n' && *p != '\r') {
            if (i < maxlen - 1) out[i++] = *p;
            p++;
        }
        if (*p == ',') p++;
    }
    out[i] = '\0';
    return p;
}

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_hyg.csv> <output.csv> [mag_limit]\n", argv[0]);
        return 1;
    }

    const char *infile = argv[1];
    const char *outfile = argv[2];
    float mag_limit = 6.5f;
    if (argc >= 4)
        mag_limit = (float)atof(argv[3]);

    FILE *fin = fopen(infile, "r");
    if (!fin) { perror("Cannot open input"); return 1; }

    FILE *fout = fopen(outfile, "w");
    if (!fout) { perror("Cannot open output"); fclose(fin); return 1; }

    char line[MAX_LINE];
    char field[256];
    int count = 0;
    int linenum = 0;

    fprintf(fout, "hip_id,ra,dec,mag\n");

    while (fgets(line, sizeof(line), fin)) {
        linenum++;
        if (linenum == 1) continue; /* skip header */

        const char *p = line;

        /* Field 0: id */
        p = next_field(p, field, sizeof(field));

        /* Field 1: hip */
        p = next_field(p, field, sizeof(field));
        int hip = atoi(field);
        if (hip <= 0) continue;

        /* Field 2: hd */
        p = next_field(p, field, sizeof(field));
        /* Field 3: hr */
        p = next_field(p, field, sizeof(field));
        /* Field 4: gl */
        p = next_field(p, field, sizeof(field));
        /* Field 5: bf */
        p = next_field(p, field, sizeof(field));
        /* Field 6: proper */
        p = next_field(p, field, sizeof(field));

        /* Field 7: ra (hours) */
        p = next_field(p, field, sizeof(field));
        float ra_hours = (float)atof(field);
        float ra_deg = ra_hours * 15.0f;  /* convert hours to degrees */

        /* Field 8: dec (degrees) */
        p = next_field(p, field, sizeof(field));
        float dec_deg = (float)atof(field);

        /* Skip fields 9-12: dist, pmra, pmdec, rv */
        p = next_field(p, field, sizeof(field));
        p = next_field(p, field, sizeof(field));
        p = next_field(p, field, sizeof(field));
        p = next_field(p, field, sizeof(field));

        /* Field 13: mag */
        p = next_field(p, field, sizeof(field));
        float mag = (float)atof(field);

        if (mag > mag_limit) continue;

        fprintf(fout, "%d,%.6f,%.6f,%.2f\n", hip, ra_deg, dec_deg, mag);
        count++;
    }

    fclose(fin);
    fclose(fout);

    printf("Converted %d Hipparcos stars (mag <= %.1f) -> %s\n", count, mag_limit, outfile);
    return 0;
}
