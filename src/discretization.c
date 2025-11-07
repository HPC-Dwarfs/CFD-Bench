/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "allocate.h"
#include "discretization.h"
#include "parameter.h"
#include "util.h"

#define P(i, j, k)                                                                       \
    p[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define F(i, j, k)                                                                       \
    f[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define G(i, j, k)                                                                       \
    g[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define H(i, j, k)                                                                       \
    h[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define U(i, j, k)                                                                       \
    u[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define V(i, j, k)                                                                       \
    v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define W(i, j, k)                                                                       \
    w[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define RHS(i, j, k)                                                                     \
    rhs[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

static void printConfig(Discretization* s)
{
    if (commIsMaster(&s->comm)) {
        printf("Parameters for #%s#\n", s->problem);
        printf("BC Left:%d Right:%d Bottom:%d Top:%d Front:%d Back:%d\n",
            s->bcLeft,
            s->bcRight,
            s->bcBottom,
            s->bcTop,
            s->bcFront,
            s->bcBack);
        printf("\tReynolds number: %.2f\n", s->re);
        printf("\tGx Gy: %.2f %.2f %.2f\n", s->gx, s->gy, s->gz);
        printf("Geometry data:\n");
        printf("\tDomain box size (x, y, z): %.2f, %.2f, %.2f\n",
            s->grid.xlength,
            s->grid.ylength,
            s->grid.zlength);
        printf("\tCells (x, y, z): %d, %d, %d\n",
            s->grid.imax,
            s->grid.jmax,
            s->grid.kmax);
        printf("\tCell size (dx, dy, dz): %f, %f, %f\n",
            s->grid.dx,
            s->grid.dy,
            s->grid.dz);
        printf("Timestep parameters:\n");
        printf("\tDefault stepsize: %.2f, Final time %.2f\n", s->dt, s->te);
        printf("\tdt bound: %.6f\n", s->dtBound);
        printf("\tTau factor: %.2f\n", s->tau);
        printf("Iterative parameters:\n");
        printf("\tMax iterations: %d\n", s->itermax);
        printf("\tepsilon (stopping tolerance) : %f\n", s->eps);
        printf("\tgamma factor: %f\n", s->gamma);
        printf("\tomega (SOR relaxation): %f\n", s->omega);
    }
    commPrintConfig(&s->comm);
}

void initDiscretization(Discretization* s, Parameter* params)
{
    s->problem  = params->name;
    s->bcLeft   = params->bcLeft;
    s->bcRight  = params->bcRight;
    s->bcBottom = params->bcBottom;
    s->bcTop    = params->bcTop;
    s->bcFront  = params->bcFront;
    s->bcBack   = params->bcBack;

    s->grid.imax    = params->imax;
    s->grid.jmax    = params->jmax;
    s->grid.kmax    = params->kmax;
    s->grid.xlength = params->xlength;
    s->grid.ylength = params->ylength;
    s->grid.zlength = params->zlength;
    s->grid.dx      = params->xlength / params->imax;
    s->grid.dy      = params->ylength / params->jmax;
    s->grid.dz      = params->zlength / params->kmax;

    s->eps     = params->eps;
    s->omega   = params->omg;
    s->itermax = params->itermax;
    s->re      = params->re;
    s->gx      = params->gx;
    s->gy      = params->gy;
    s->gz      = params->gz;
    s->dt      = params->dt;
    s->te      = params->te;
    s->tau     = params->tau;
    s->gamma   = params->gamma;

    /* allocate arrays */
    int imaxLocal = s->comm.imaxLocal;
    int jmaxLocal = s->comm.jmaxLocal;
    int kmaxLocal = s->comm.kmaxLocal;
    size_t size   = (imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

    s->u   = allocate(64, size * sizeof(double));
    s->v   = allocate(64, size * sizeof(double));
    s->w   = allocate(64, size * sizeof(double));
    s->p   = allocate(64, size * sizeof(double));
    s->rhs = allocate(64, size * sizeof(double));
    s->f   = allocate(64, size * sizeof(double));
    s->g   = allocate(64, size * sizeof(double));
    s->h   = allocate(64, size * sizeof(double));

    for (int i = 0; i < size; i++) {
        s->u[i]   = params->u_init;
        s->v[i]   = params->v_init;
        s->w[i]   = params->w_init;
        s->p[i]   = params->p_init;
        s->rhs[i] = 0.0;
        s->f[i]   = 0.0;
        s->g[i]   = 0.0;
        s->h[i]   = 0.0;
    }

    double dx = s->grid.dx;
    double dy = s->grid.dy;
    double dz = s->grid.dz;

    double invSqrSum = 1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
    s->dtBound       = 0.5 * s->re * 1.0 / invSqrSum;

#ifdef VERBOSE
    printConfig(s);
#endif /* VERBOSE */
}

void setBoundaryConditions(Discretization* s)
{
    int imaxLocal = s->comm.imaxLocal;
    int jmaxLocal = s->comm.jmaxLocal;
    int kmaxLocal = s->comm.kmaxLocal;

    double* u = s->u;
    double* v = s->v;
    double* w = s->w;

    if (commIsBoundary(&s->comm, TOP)) {
        switch (s->bcTop) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, jmaxLocal + 1, k) = -U(i, jmaxLocal, k);
                    V(i, jmaxLocal, k)     = 0.0;
                    W(i, jmaxLocal + 1, k) = -W(i, jmaxLocal, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, jmaxLocal + 1, k) = U(i, jmaxLocal, k);
                    V(i, jmaxLocal, k)     = 0.0;
                    W(i, jmaxLocal + 1, k) = W(i, jmaxLocal, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, jmaxLocal + 1, k) = U(i, jmaxLocal, k);
                    V(i, jmaxLocal, k)     = V(i, jmaxLocal - 1, k);
                    W(i, jmaxLocal + 1, k) = W(i, jmaxLocal, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&s->comm, BOTTOM)) {
        switch (s->bcBottom) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, 0, k) = -U(i, 1, k);
                    V(i, 0, k) = 0.0;
                    W(i, 0, k) = -W(i, 1, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, 0, k) = U(i, 1, k);
                    V(i, 0, k) = 0.0;
                    W(i, 0, k) = W(i, 1, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, 0, k) = U(i, 1, k);
                    V(i, 0, k) = V(i, 1, k);
                    W(i, 0, k) = W(i, 1, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&s->comm, LEFT)) {
        switch (s->bcLeft) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = 0.0;
                    V(0, j, k) = -V(1, j, k);
                    W(0, j, k) = -W(1, j, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = 0.0;
                    V(0, j, k) = V(1, j, k);
                    W(0, j, k) = W(1, j, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = U(1, j, k);
                    V(0, j, k) = V(1, j, k);
                    W(0, j, k) = W(1, j, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&s->comm, RIGHT)) {
        switch (s->bcRight) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(imaxLocal, j, k)     = 0.0;
                    V(imaxLocal + 1, j, k) = -V(imaxLocal, j, k);
                    W(imaxLocal + 1, j, k) = -W(imaxLocal, j, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(imaxLocal, j, k)     = 0.0;
                    V(imaxLocal + 1, j, k) = V(imaxLocal, j, k);
                    W(imaxLocal + 1, j, k) = W(imaxLocal, j, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(imaxLocal, j, k)     = U(imaxLocal - 1, j, k);
                    V(imaxLocal + 1, j, k) = V(imaxLocal, j, k);
                    W(imaxLocal + 1, j, k) = W(imaxLocal, j, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&s->comm, FRONT)) {
        switch (s->bcFront) {
        case NOSLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, 0) = -U(i, j, 1);
                    V(i, j, 0) = -V(i, j, 1);
                    W(i, j, 0) = 0.0;
                }
            }
            break;
        case SLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, 0) = U(i, j, 1);
                    V(i, j, 0) = V(i, j, 1);
                    W(i, j, 0) = 0.0;
                }
            }
            break;
        case OUTFLOW:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, 0) = U(i, j, 1);
                    V(i, j, 0) = V(i, j, 1);
                    W(i, j, 0) = W(i, j, 1);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&s->comm, BACK)) {
        switch (s->bcBack) {
        case NOSLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, kmaxLocal + 1) = -U(i, j, kmaxLocal);
                    V(i, j, kmaxLocal + 1) = -V(i, j, kmaxLocal);
                    W(i, j, kmaxLocal)     = 0.0;
                }
            }
            break;
        case SLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, kmaxLocal + 1) = U(i, j, kmaxLocal);
                    V(i, j, kmaxLocal + 1) = V(i, j, kmaxLocal);
                    W(i, j, kmaxLocal)     = 0.0;
                }
            }
            break;
        case OUTFLOW:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, kmaxLocal + 1) = U(i, j, kmaxLocal);
                    V(i, j, kmaxLocal + 1) = V(i, j, kmaxLocal);
                    W(i, j, kmaxLocal)     = W(i, j, kmaxLocal - 1);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }
}

void computeRHS(Discretization* s)
{
    int imaxLocal = s->comm.imaxLocal;
    int jmaxLocal = s->comm.jmaxLocal;
    int kmaxLocal = s->comm.kmaxLocal;

    double idx = 1.0 / s->grid.dx;
    double idy = 1.0 / s->grid.dy;
    double idz = 1.0 / s->grid.dz;
    double idt = 1.0 / s->dt;

    double* rhs = s->rhs;
    double* f   = s->f;
    double* g   = s->g;
    double* h   = s->h;

    commShift(&s->comm, f, g, h);

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                RHS(i, j, k) = ((F(i, j, k) - F(i - 1, j, k)) * idx +
                                   (G(i, j, k) - G(i, j - 1, k)) * idy +
                                   (H(i, j, k) - H(i, j, k - 1)) * idz) *
                               idt;
            }
        }
    }
}

void setSpecialBoundaryCondition(Discretization* s)
{
    int imaxLocal = s->comm.imaxLocal;
    int jmaxLocal = s->comm.jmaxLocal;
    int kmaxLocal = s->comm.kmaxLocal;

    double* u = s->u;

    if (strcmp(s->problem, "dcavity") == 0) {
        if (commIsBoundary(&s->comm, TOP)) {
            for (int k = 1; k < kmaxLocal; k++) {
                for (int i = 1; i < imaxLocal; i++) {
                    U(i, jmaxLocal + 1, k) = 2.0 - U(i, jmaxLocal, k);
                }
            }
        }
    } else if (strcmp(s->problem, "canal") == 0) {
        if (commIsBoundary(&s->comm, LEFT)) {
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = 2.0;
                }
            }
        }
    }
}

static double maxElement(Discretization* s, double* m)
{
    int size = (s->comm.imaxLocal + 2) * (s->comm.jmaxLocal + 2) *
               (s->comm.kmaxLocal + 2);
    double maxval = DBL_MIN;

    for (int i = 0; i < size; i++) {
        maxval = MAX(maxval, fabs(m[i]));
    }
    commReduction(&maxval, MAX);
    return maxval;
}

void normalizePressure(Discretization* s)
{
    int imaxLocal = s->comm.imaxLocal;
    int jmaxLocal = s->comm.jmaxLocal;
    int kmaxLocal = s->comm.kmaxLocal;

    double* p   = s->p;
    double avgP = 0.0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                avgP += P(i, j, k);
            }
        }
    }
    commReduction(&avgP, SUM);
    avgP /= (s->grid.imax * s->grid.jmax * s->grid.kmax);

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                P(i, j, k) = P(i, j, k) - avgP;
            }
        }
    }
}

void computeTimestep(Discretization* s)
{
    double dt = s->dtBound;
    double dx = s->grid.dx;
    double dy = s->grid.dy;
    double dz = s->grid.dz;

    double umax = maxElement(s, s->u);
    double vmax = maxElement(s, s->v);
    double wmax = maxElement(s, s->w);

    if (umax > 0) {
        dt = (dt > dx / umax) ? dx / umax : dt;
    }
    if (vmax > 0) {
        dt = (dt > dy / vmax) ? dy / vmax : dt;
    }
    if (wmax > 0) {
        dt = (dt > dz / wmax) ? dz / wmax : dt;
    }

    s->dt = dt * s->tau;
}



void computeFG(Discretization* s)
{
    int imaxLocal = s->comm.imaxLocal;
    int jmaxLocal = s->comm.jmaxLocal;
    int kmaxLocal = s->comm.kmaxLocal;

    double* u = s->u;
    double* v = s->v;
    double* w = s->w;
    double* f = s->f;
    double* g = s->g;
    double* h = s->h;

    double gx = s->gx;
    double gy = s->gy;
    double gz = s->gz;
    double dt = s->dt;

    double gamma     = s->gamma;
    double inverseRe = 1.0 / s->re;
    double inverseDx = 1.0 / s->grid.dx;
    double inverseDy = 1.0 / s->grid.dy;
    double inverseDz = 1.0 / s->grid.dz;
    double du2dx, dv2dy, dw2dz;
    double duvdx, duwdx, duvdy, dvwdy, duwdz, dvwdz;
    double du2dx2, du2dy2, du2dz2;
    double dv2dx2, dv2dy2, dv2dz2;
    double dw2dx2, dw2dy2, dw2dz2;

    commExchange(&s->comm, u);
    commExchange(&s->comm, v);
    commExchange(&s->comm, w);

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                du2dx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i + 1, j, k)) *
                                    (U(i, j, k) + U(i + 1, j, k)) -
                                (U(i, j, k) + U(i - 1, j, k)) *
                                    (U(i, j, k) + U(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i + 1, j, k)) *
                                    (U(i, j, k) - U(i + 1, j, k)) +
                                fabs(U(i, j, k) + U(i - 1, j, k)) *
                                    (U(i, j, k) - U(i - 1, j, k)));

                duvdy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i + 1, j, k)) *
                                    (U(i, j, k) + U(i, j + 1, k)) -
                                (V(i, j - 1, k) + V(i + 1, j - 1, k)) *
                                    (U(i, j, k) + U(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i + 1, j, k)) *
                                    (U(i, j, k) - U(i, j + 1, k)) +
                                fabs(V(i, j - 1, k) + V(i + 1, j - 1, k)) *
                                    (U(i, j, k) - U(i, j - 1, k)));

                duwdz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i + 1, j, k)) *
                                    (U(i, j, k) + U(i, j, k + 1)) -
                                (W(i, j, k - 1) + W(i + 1, j, k - 1)) *
                                    (U(i, j, k) + U(i, j, k - 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i + 1, j, k)) *
                                    (U(i, j, k) - U(i, j, k + 1)) +
                                fabs(W(i, j, k - 1) + W(i + 1, j, k - 1)) *
                                    (U(i, j, k) - U(i, j, k - 1)));

                du2dx2 = inverseDx * inverseDx *
                         (U(i + 1, j, k) - 2.0 * U(i, j, k) + U(i - 1, j, k));
                du2dy2 = inverseDy * inverseDy *
                         (U(i, j + 1, k) - 2.0 * U(i, j, k) + U(i, j - 1, k));
                du2dz2 = inverseDz * inverseDz *
                         (U(i, j, k + 1) - 2.0 * U(i, j, k) + U(i, j, k - 1));
                F(i, j, k) = U(i, j, k) + dt * (inverseRe * (du2dx2 + du2dy2 + du2dz2) -
                                                   du2dx - duvdy - duwdz + gx);

                duvdx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i, j + 1, k)) *
                                    (V(i, j, k) + V(i + 1, j, k)) -
                                (U(i - 1, j, k) + U(i - 1, j + 1, k)) *
                                    (V(i, j, k) + V(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i, j + 1, k)) *
                                    (V(i, j, k) - V(i + 1, j, k)) +
                                fabs(U(i - 1, j, k) + U(i - 1, j + 1, k)) *
                                    (V(i, j, k) - V(i - 1, j, k)));

                dv2dy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i, j + 1, k)) *
                                    (V(i, j, k) + V(i, j + 1, k)) -
                                (V(i, j, k) + V(i, j - 1, k)) *
                                    (V(i, j, k) + V(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i, j + 1, k)) *
                                    (V(i, j, k) - V(i, j + 1, k)) +
                                fabs(V(i, j, k) + V(i, j - 1, k)) *
                                    (V(i, j, k) - V(i, j - 1, k)));

                dvwdz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i, j + 1, k)) *
                                    (V(i, j, k) + V(i, j, k + 1)) -
                                (W(i, j, k - 1) + W(i, j + 1, k - 1)) *
                                    (V(i, j, k) + V(i, j, k + 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i, j + 1, k)) *
                                    (V(i, j, k) - V(i, j, k + 1)) +
                                fabs(W(i, j, k - 1) + W(i, j + 1, k - 1)) *
                                    (V(i, j, k) - V(i, j, k + 1)));

                dv2dx2 = inverseDx * inverseDx *
                         (V(i + 1, j, k) - 2.0 * V(i, j, k) + V(i - 1, j, k));
                dv2dy2 = inverseDy * inverseDy *
                         (V(i, j + 1, k) - 2.0 * V(i, j, k) + V(i, j - 1, k));
                dv2dz2 = inverseDz * inverseDz *
                         (V(i, j, k + 1) - 2.0 * V(i, j, k) + V(i, j, k - 1));
                G(i, j, k) = V(i, j, k) + dt * (inverseRe * (dv2dx2 + dv2dy2 + dv2dz2) -
                                                   duvdx - dv2dy - dvwdz + gy);

                duwdx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i, j, k + 1)) *
                                    (W(i, j, k) + W(i + 1, j, k)) -
                                (U(i - 1, j, k) + U(i - 1, j, k + 1)) *
                                    (W(i, j, k) + W(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i, j, k + 1)) *
                                    (W(i, j, k) - W(i + 1, j, k)) +
                                fabs(U(i - 1, j, k) + U(i - 1, j, k + 1)) *
                                    (W(i, j, k) - W(i - 1, j, k)));

                dvwdy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i, j, k + 1)) *
                                    (W(i, j, k) + W(i, j + 1, k)) -
                                (V(i, j - 1, k + 1) + V(i, j - 1, k)) *
                                    (W(i, j, k) + W(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i, j, k + 1)) *
                                    (W(i, j, k) - W(i, j + 1, k)) +
                                fabs(V(i, j - 1, k + 1) + V(i, j - 1, k)) *
                                    (W(i, j, k) - W(i, j - 1, k)));

                dw2dz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i, j, k + 1)) *
                                    (W(i, j, k) + W(i, j, k + 1)) -
                                (W(i, j, k) + W(i, j, k - 1)) *
                                    (W(i, j, k) + W(i, j, k - 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i, j, k + 1)) *
                                    (W(i, j, k) - W(i, j, k + 1)) +
                                fabs(W(i, j, k) + W(i, j, k - 1)) *
                                    (W(i, j, k) - W(i, j, k - 1)));

                dw2dx2 = inverseDx * inverseDx *
                         (W(i + 1, j, k) - 2.0 * W(i, j, k) + W(i - 1, j, k));
                dw2dy2 = inverseDy * inverseDy *
                         (W(i, j + 1, k) - 2.0 * W(i, j, k) + W(i, j - 1, k));
                dw2dz2 = inverseDz * inverseDz *
                         (W(i, j, k + 1) - 2.0 * W(i, j, k) + W(i, j, k - 1));
                H(i, j, k) = W(i, j, k) + dt * (inverseRe * (dw2dx2 + dw2dy2 + dw2dz2) -
                                                   duwdx - dvwdy - dw2dz + gz);
            }
        }
    }

    /* ----------------------------- boundary of F ---------------------------
     */
    if (commIsBoundary(&s->comm, LEFT)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int j = 1; j < jmaxLocal + 1; j++) {
                F(0, j, k) = U(0, j, k);
            }
        }
    }

    if (commIsBoundary(&s->comm, RIGHT)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int j = 1; j < jmaxLocal + 1; j++) {
                F(imaxLocal, j, k) = U(imaxLocal, j, k);
            }
        }
    }

    /* ----------------------------- boundary of G ---------------------------
     */
    if (commIsBoundary(&s->comm, BOTTOM)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                G(i, 0, k) = V(i, 0, k);
            }
        }
    }

    if (commIsBoundary(&s->comm, TOP)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                G(i, jmaxLocal, k) = V(i, jmaxLocal, k);
            }
        }
    }

    /* ----------------------------- boundary of H ---------------------------
     */
    if (commIsBoundary(&s->comm, FRONT)) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                H(i, j, 0) = W(i, j, 0);
            }
        }
    }

    if (commIsBoundary(&s->comm, BACK)) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                H(i, j, kmaxLocal) = W(i, j, kmaxLocal);
            }
        }
    }
}

void adaptUV(Discretization* s)
{
    int imaxLocal = s->comm.imaxLocal;
    int jmaxLocal = s->comm.jmaxLocal;
    int kmaxLocal = s->comm.kmaxLocal;

    double* p = s->p;
    double* u = s->u;
    double* v = s->v;
    double* w = s->w;
    double* f = s->f;
    double* g = s->g;
    double* h = s->h;

    double factorX = s->dt / s->grid.dx;
    double factorY = s->dt / s->grid.dy;
    double factorZ = s->dt / s->grid.dz;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                U(i, j, k) = F(i, j, k) - (P(i + 1, j, k) - P(i, j, k)) * factorX;
                V(i, j, k) = G(i, j, k) - (P(i, j + 1, k) - P(i, j, k)) * factorY;
                W(i, j, k) = H(i, j, k) - (P(i, j, k + 1) - P(i, j, k)) * factorZ;
            }
        }
    }
}