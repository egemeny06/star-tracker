#include "attitude.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define RAD2DEG(r) ((r) * (float)(180.0 / M_PI))

/* ------------------------------------------------------------------ */
/*  Helper: cross product  c = a x b                                  */
/* ------------------------------------------------------------------ */
static void cross(const Vec3 *a, const Vec3 *b, Vec3 *c)
{
    c->x = a->y * b->z - a->z * b->y;
    c->y = a->z * b->x - a->x * b->z;
    c->z = a->x * b->y - a->y * b->x;
}

/* ------------------------------------------------------------------ */
/*  Helper: normalize vector in-place                                 */
/* ------------------------------------------------------------------ */
static void vec3_normalize(Vec3 *v)
{
    float n = sqrtf(v->x * v->x + v->y * v->y + v->z * v->z);
    if (n > 1e-12f) {
        v->x /= n;
        v->y /= n;
        v->z /= n;
    }
}

/* ================================================================== */
/*  Match list                                                        */
/* ================================================================== */

void match_list_init(MatchList *ml)
{
    memset(ml, 0, sizeof(MatchList));
}

int match_list_add(MatchList *ml,
                   float bx, float by, float bz,
                   float rx, float ry, float rz,
                   float weight)
{
    if (ml->count >= MAX_MATCHED_STARS)
        return -1;

    StarMatch *sm = &ml->matches[ml->count];
    sm->body.x = bx;  sm->body.y = by;  sm->body.z = bz;
    sm->ref.x  = rx;  sm->ref.y  = ry;  sm->ref.z  = rz;
    sm->weight = weight;
    ml->count++;
    return 0;
}

/* ================================================================== */
/*  TRIAD Algorithm                                                   */
/*  Builds two triads (orthonormal bases) from body and reference     */
/*  vectors, then computes DCM = M_body * M_ref^T                    */
/* ================================================================== */

int attitude_triad(const MatchList *ml, DCM *dcm)
{
    if (ml->count < 2)
        return -1;

    Vec3 b1 = ml->matches[0].body;
    Vec3 b2 = ml->matches[1].body;
    Vec3 r1 = ml->matches[0].ref;
    Vec3 r2 = ml->matches[1].ref;

    vec3_normalize(&b1);
    vec3_normalize(&b2);
    vec3_normalize(&r1);
    vec3_normalize(&r2);

    /* Body triad */
    Vec3 tb1 = b1;                     /* t1 = b1 */
    Vec3 tb2;
    cross(&b1, &b2, &tb2);             /* t2 = b1 x b2 */
    vec3_normalize(&tb2);
    Vec3 tb3;
    cross(&tb1, &tb2, &tb3);           /* t3 = t1 x t2 */

    /* Reference triad */
    Vec3 tr1 = r1;                     /* s1 = r1 */
    Vec3 tr2;
    cross(&r1, &r2, &tr2);             /* s2 = r1 x r2 */
    vec3_normalize(&tr2);
    Vec3 tr3;
    cross(&tr1, &tr2, &tr3);           /* s3 = s1 x s2 */

    /* DCM = [tb1 tb2 tb3] * [tr1 tr2 tr3]^T
       A_ij = sum_k M_body(i,k) * M_ref(j,k)   */
    float Mb[3][3] = {
        {tb1.x, tb2.x, tb3.x},
        {tb1.y, tb2.y, tb3.y},
        {tb1.z, tb2.z, tb3.z}
    };
    float Mr[3][3] = {
        {tr1.x, tr2.x, tr3.x},
        {tr1.y, tr2.y, tr3.y},
        {tr1.z, tr2.z, tr3.z}
    };

    /* C = Mb * Mr^T */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dcm->m[i][j] = 0.0f;
            for (int k = 0; k < 3; k++) {
                dcm->m[i][j] += Mb[i][k] * Mr[j][k];
            }
        }
    }

    return 0;
}

/* ================================================================== */
/*  QUEST Algorithm                                                   */
/*  Davenport's q-method solved via Newton-Raphson on characteristic  */
/*  equation. Embedded-friendly: no eigenvalue decomposition needed.  */
/* ================================================================== */

int attitude_quest(const MatchList *ml, Quaternion *q)
{
    if (ml->count < 2)
        return -1;

    uint16_t n = ml->count;

    /* --- Step 1: Build B matrix = sum( w_i * b_i * r_i^T ) --- */
    float B[3][3];
    memset(B, 0, sizeof(B));

    float sum_w = 0.0f;
    for (uint16_t k = 0; k < n; k++) {
        float w = ml->matches[k].weight;
        const Vec3 *b = &ml->matches[k].body;
        const Vec3 *r = &ml->matches[k].ref;
        sum_w += w;

        B[0][0] += w * b->x * r->x;
        B[0][1] += w * b->x * r->y;
        B[0][2] += w * b->x * r->z;
        B[1][0] += w * b->y * r->x;
        B[1][1] += w * b->y * r->y;
        B[1][2] += w * b->y * r->z;
        B[2][0] += w * b->z * r->x;
        B[2][1] += w * b->z * r->y;
        B[2][2] += w * b->z * r->z;
    }

    /* --- Step 2: Compute S = B + B^T, sigma = trace(B) --- */
    float S[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            S[i][j] = B[i][j] + B[j][i];

    float sigma = B[0][0] + B[1][1] + B[2][2];

    /* --- Step 3: Z vector (antisymmetric part of B) ---
       Convention: B = sum(w * b * r^T) gives C_br (body from ref).
       Z must match: Z_i = (B^T - B) antisymmetric extraction */
    Vec3 Z;
    Z.x = B[2][1] - B[1][2];
    Z.y = B[0][2] - B[2][0];
    Z.z = B[1][0] - B[0][1];

    /* --- Step 4: Newton-Raphson to find max eigenvalue lambda --- */
    /* Characteristic equation of K matrix (4x4 Davenport):
       f(l) = l^4 - (a+b)*l^2 - c*l + (a*b + c*sigma - delta) = 0
       where a = sigma^2 - kappa, b = sigma^2 + |Z|^2
       Simplified: f(l) = (l^2 - a)(l - sigma) - delta - |Z|^2 * (l + sigma)  ... nope
       Use direct form:
       f(l) = l^4 - (a+b)*l^2 - c*l + (a*b + c*sigma - delta)

       Actually the standard QUEST characteristic equation is:
       f(l) = (l - sigma)*(l + sigma)*(l^2 - a) - delta*(l - sigma)
            - Z^T * S * Z * (l + sigma) ... this gets complicated.

       Let's use the simpler Shuster form directly:
       det[(l+sigma)*I - S] * (l - sigma) - Z^T * adj[(l+sigma)*I - S] * Z = 0
    */

    /* Start from lambda_0 = sum_w (upper bound for max eigenvalue) */
    float lambda = sum_w;

    for (int iter = 0; iter < 100; iter++) {
        /* Matrix A = (lambda + sigma)*I - S */
        float A00 = lambda + sigma - S[0][0];
        float A11 = lambda + sigma - S[1][1];
        float A22 = lambda + sigma - S[2][2];
        float A01 = -S[0][1];
        float A02 = -S[0][2];
        float A12 = -S[1][2];

        /* Cofactors of A (adjugate) for symmetric matrix */
        float adj00 = A11 * A22 - A12 * A12;
        float adj01 = A02 * A12 - A01 * A22;
        float adj02 = A01 * A12 - A02 * A11;
        float adj11 = A00 * A22 - A02 * A02;
        float adj12 = A02 * A01 - A00 * A12;
        float adj22 = A00 * A11 - A01 * A01;

        /* det(A) */
        float detA = A00 * adj00 + A01 * adj01 + A02 * adj02;

        /* Z^T * adj(A) * Z */
        float ZtAdjZ = Z.x * (adj00*Z.x + adj01*Z.y + adj02*Z.z)
                     + Z.y * (adj01*Z.x + adj11*Z.y + adj12*Z.z)
                     + Z.z * (adj02*Z.x + adj12*Z.y + adj22*Z.z);

        /* f(lambda) = detA * (lambda - sigma) - ZtAdjZ */
        float f_val = detA * (lambda - sigma) - ZtAdjZ;

        /* f'(lambda) via numerical derivative */
        float eps = 1e-6f;
        float lp = lambda + eps;
        float A00p = lp + sigma - S[0][0];
        float A11p = lp + sigma - S[1][1];
        float A22p = lp + sigma - S[2][2];

        float adj00p = A11p * A22p - A12 * A12;
        float adj01p = A02 * A12 - A01 * A22p;
        float adj02p = A01 * A12 - A02 * A11p;
        float adj11p = A00p * A22p - A02 * A02;
        float adj12p = A02 * A01 - A00p * A12;
        float adj22p = A00p * A11p - A01 * A01;

        float detAp = A00p * adj00p + A01 * adj01p + A02 * adj02p;

        float ZtAdjZp = Z.x * (adj00p*Z.x + adj01p*Z.y + adj02p*Z.z)
                      + Z.y * (adj01p*Z.x + adj11p*Z.y + adj12p*Z.z)
                      + Z.z * (adj02p*Z.x + adj12p*Z.y + adj22p*Z.z);

        float f_valp = detAp * (lp - sigma) - ZtAdjZp;
        float df = (f_valp - f_val) / eps;

        if (fabsf(df) < 1e-15f)
            break;

        float lambda_new = lambda - f_val / df;
        if (fabsf(lambda_new - lambda) < 1e-10f) {
            lambda = lambda_new;
            break;
        }
        lambda = lambda_new;
    }

    /* --- Step 6: Compute quaternion from converged lambda --- */
    {
        float A00 = lambda + sigma - S[0][0];
        float A11 = lambda + sigma - S[1][1];
        float A22 = lambda + sigma - S[2][2];
        float A01 = -S[0][1];
        float A02 = -S[0][2];
        float A12 = -S[1][2];

        /* Adjugate of A */
        float adj00 = A11 * A22 - A12 * A12;
        float adj01 = A02 * A12 - A01 * A22;
        float adj02 = A01 * A12 - A02 * A11;
        float adj11 = A00 * A22 - A02 * A02;
        float adj12 = A02 * A01 - A00 * A12;
        float adj22 = A00 * A11 - A01 * A01;

        /* Rodrigues vector: p = adj(A) * Z / det(A)
           But for QUEST we use: q_vec = adj(A) * Z, q_scalar = det(A) */
        float qx = adj00 * Z.x + adj01 * Z.y + adj02 * Z.z;
        float qy = adj01 * Z.x + adj11 * Z.y + adj12 * Z.z;
        float qz = adj02 * Z.x + adj12 * Z.y + adj22 * Z.z;
        float detA = A00 * adj00 + A01 * adj01 + A02 * adj02;

        /* Quaternion (unnormalized): [qx, qy, qz, detA] */
        q->q1 = qx;
        q->q2 = qy;
        q->q3 = qz;
        q->q4 = detA;
    }

    quaternion_normalize(q);
    return 0;
}

/* ================================================================== */
/*  Conversions                                                       */
/* ================================================================== */

void quaternion_normalize(Quaternion *q)
{
    float n = sqrtf(q->q1*q->q1 + q->q2*q->q2 + q->q3*q->q3 + q->q4*q->q4);
    if (n > 1e-12f) {
        q->q1 /= n;
        q->q2 /= n;
        q->q3 /= n;
        q->q4 /= n;
    }
    /* Enforce positive scalar part for unique representation */
    if (q->q4 < 0.0f) {
        q->q1 = -q->q1;
        q->q2 = -q->q2;
        q->q3 = -q->q3;
        q->q4 = -q->q4;
    }
}

void dcm_to_quaternion(const DCM *dcm, Quaternion *q)
{
    const float (*m)[3] = dcm->m;
    float tr = m[0][0] + m[1][1] + m[2][2];

    if (tr > 0.0f) {
        float s = 2.0f * sqrtf(tr + 1.0f);
        q->q4 = 0.25f * s;
        q->q1 = (m[2][1] - m[1][2]) / s;
        q->q2 = (m[0][2] - m[2][0]) / s;
        q->q3 = (m[1][0] - m[0][1]) / s;
    } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
        float s = 2.0f * sqrtf(1.0f + m[0][0] - m[1][1] - m[2][2]);
        q->q4 = (m[2][1] - m[1][2]) / s;
        q->q1 = 0.25f * s;
        q->q2 = (m[0][1] + m[1][0]) / s;
        q->q3 = (m[0][2] + m[2][0]) / s;
    } else if (m[1][1] > m[2][2]) {
        float s = 2.0f * sqrtf(1.0f + m[1][1] - m[0][0] - m[2][2]);
        q->q4 = (m[0][2] - m[2][0]) / s;
        q->q1 = (m[0][1] + m[1][0]) / s;
        q->q2 = 0.25f * s;
        q->q3 = (m[1][2] + m[2][1]) / s;
    } else {
        float s = 2.0f * sqrtf(1.0f + m[2][2] - m[0][0] - m[1][1]);
        q->q4 = (m[1][0] - m[0][1]) / s;
        q->q1 = (m[0][2] + m[2][0]) / s;
        q->q2 = (m[1][2] + m[2][1]) / s;
        q->q3 = 0.25f * s;
    }
    quaternion_normalize(q);
}

void quaternion_to_dcm(const Quaternion *q, DCM *dcm)
{
    float q1 = q->q1, q2 = q->q2, q3 = q->q3, q4 = q->q4;

    dcm->m[0][0] = 1.0f - 2.0f * (q2*q2 + q3*q3);
    dcm->m[0][1] = 2.0f * (q1*q2 - q3*q4);
    dcm->m[0][2] = 2.0f * (q1*q3 + q2*q4);

    dcm->m[1][0] = 2.0f * (q1*q2 + q3*q4);
    dcm->m[1][1] = 1.0f - 2.0f * (q1*q1 + q3*q3);
    dcm->m[1][2] = 2.0f * (q2*q3 - q1*q4);

    dcm->m[2][0] = 2.0f * (q1*q3 - q2*q4);
    dcm->m[2][1] = 2.0f * (q2*q3 + q1*q4);
    dcm->m[2][2] = 1.0f - 2.0f * (q1*q1 + q2*q2);
}

void dcm_to_euler(const DCM *dcm, float *roll, float *pitch, float *yaw)
{
    const float (*m)[3] = dcm->m;

    *pitch = RAD2DEG(asinf(-m[2][0]));

    if (fabsf(m[2][0]) < 0.9999f) {
        *roll = RAD2DEG(atan2f(m[2][1], m[2][2]));
        *yaw  = RAD2DEG(atan2f(m[1][0], m[0][0]));
    } else {
        /* Gimbal lock */
        *roll = RAD2DEG(atan2f(-m[0][1], m[1][1]));
        *yaw  = 0.0f;
    }
}

/* ================================================================== */
/*  Utilities                                                         */
/* ================================================================== */

void dcm_multiply(const DCM *a, const DCM *b, DCM *c)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            c->m[i][j] = 0.0f;
            for (int k = 0; k < 3; k++)
                c->m[i][j] += a->m[i][k] * b->m[k][j];
        }
}

void dcm_transpose(const DCM *in, DCM *out)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            out->m[i][j] = in->m[j][i];
}

float quaternion_error_angle(const Quaternion *q1, const Quaternion *q2)
{
    /* Error quaternion: q_err = q1^-1 * q2 */
    /* For unit quaternions, |dot(q1,q2)| gives cos(theta/2) */
    float d = q1->q1 * q2->q1 + q1->q2 * q2->q2
            + q1->q3 * q2->q3 + q1->q4 * q2->q4;

    if (d > 1.0f) d = 1.0f;
    if (d < -1.0f) d = -1.0f;

    /* Total rotation angle between two attitudes */
    return 2.0f * RAD2DEG(acosf(fabsf(d)));
}

void dcm_print(const DCM *dcm)
{
    printf("DCM (Direction Cosine Matrix):\n");
    for (int i = 0; i < 3; i++) {
        printf("  [%+.6f  %+.6f  %+.6f]\n",
               dcm->m[i][0], dcm->m[i][1], dcm->m[i][2]);
    }
}

void quaternion_print(const Quaternion *q)
{
    printf("Quaternion: [%+.6f, %+.6f, %+.6f, %+.6f] (x,y,z,w)\n",
           q->q1, q->q2, q->q3, q->q4);
}
