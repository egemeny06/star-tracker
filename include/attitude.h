#ifndef ATTITUDE_H
#define ATTITUDE_H

#include <stdint.h>

#define MAX_MATCHED_STARS 64

/* 3x3 rotation matrix (Direction Cosine Matrix) */
typedef struct {
    float m[3][3];
} DCM;

/* Quaternion: q = q4 + q1*i + q2*j + q3*k  (scalar-last convention) */
typedef struct {
    float q1;   /* vector part x */
    float q2;   /* vector part y */
    float q3;   /* vector part z */
    float q4;   /* scalar part   */
} Quaternion;

/* 3D unit vector */
typedef struct {
    float x, y, z;
} Vec3;

/* A matched star pair: body-frame observation + catalog reference */
typedef struct {
    Vec3 body;      /* unit vector in camera/body frame */
    Vec3 ref;       /* unit vector in inertial/catalog frame */
    float weight;   /* measurement weight (inverse of error) */
} StarMatch;

typedef struct {
    StarMatch matches[MAX_MATCHED_STARS];
    uint16_t count;
} MatchList;

/* ---------- Match list helpers ---------- */

void match_list_init(MatchList *ml);

int match_list_add(MatchList *ml,
                   float bx, float by, float bz,
                   float rx, float ry, float rz,
                   float weight);

/* ---------- TRIAD algorithm ---------- */
/* Uses first two matched pairs to compute attitude.
   Returns 0 on success, -1 on error. */
int attitude_triad(const MatchList *ml, DCM *dcm);

/* ---------- QUEST algorithm ---------- */
/* Optimal quaternion from N>=2 matched pairs.
   Returns 0 on success, -1 on error. */
int attitude_quest(const MatchList *ml, Quaternion *q);

/* ---------- Conversions ---------- */

void dcm_to_quaternion(const DCM *dcm, Quaternion *q);
void quaternion_to_dcm(const Quaternion *q, DCM *dcm);
void quaternion_normalize(Quaternion *q);

/* Extract Euler angles (roll, pitch, yaw) in degrees from DCM */
void dcm_to_euler(const DCM *dcm, float *roll, float *pitch, float *yaw);

/* ---------- Utilities ---------- */

void dcm_print(const DCM *dcm);
void quaternion_print(const Quaternion *q);

/* Multiply two DCMs: C = A * B */
void dcm_multiply(const DCM *a, const DCM *b, DCM *c);

/* Transpose DCM */
void dcm_transpose(const DCM *in, DCM *out);

/* Attitude error between two quaternions (angle in degrees) */
float quaternion_error_angle(const Quaternion *q1, const Quaternion *q2);

#endif /* ATTITUDE_H */
