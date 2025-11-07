/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <string.h>

#include "test.h"

#define G(v, i, j, k)                                                                    \
    v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

void testInit(Solver* s)
{
    int imaxLocal = s->comm->imaxLocal;
    int jmaxLocal = s->comm->jmaxLocal;
    int kmaxLocal = s->comm->kmaxLocal;
    int myrank    = s->comm->rank;
    double* p     = s->p;
    double* f     = s->f;
    double* g     = s->g;
    double* h     = s->h;

    for (int k = 0; k < kmaxLocal + 2; k++) {
        for (int j = 0; j < jmaxLocal + 2; j++) {
            for (int i = 0; i < imaxLocal + 2; i++) {
                G(p, i, j, k) = 10.0;
                G(f, i, j, k) = myrank + 1.0;
                G(g, i, j, k) = myrank + 1.0;
                G(h, i, j, k) = myrank + 1.0;
            }
        }
    }

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                G(p, i, j, k) = myrank + 1.0;
                G(f, i, j, k) = myrank + 1.0;
                G(g, i, j, k) = myrank + 1.0;
                G(h, i, j, k) = myrank + 1.0;
            }
        }
    }
}

static char* direction2String(Direction dir)
{
    switch (dir) {
    case LEFT:
        return "left";
        break;
    case RIGHT:
        return "right";
        break;
    case BOTTOM:
        return "bottom";
        break;
    case TOP:
        return "top";
        break;
    case FRONT:
        return "front";
        break;
    case BACK:
        return "back";
        break;
    case NDIRS:
        return "ERROR";
    default:
        return "ERROR";
    }
}

static void printPlane(Solver* s, double* a, int ymax, int xmax, Direction dir)
{
    int imaxLocal = s->comm->imaxLocal;
    int jmaxLocal = s->comm->jmaxLocal;
    int kmaxLocal = s->comm->kmaxLocal;
    char filename[50];
    snprintf(filename, 50, "halo-%s-r%d.txt", direction2String(dir), s->comm->rank);
    FILE* fh = fopen(filename, "w");

    for (int y = 0; y < ymax; y++) {
        for (int x = 0; x < xmax; x++) {
            switch (dir) {
            case LEFT:
                fprintf(fh, "%12.8f ", G(a, 0, x, y));
                break;
            case RIGHT:
                fprintf(fh, "%12.8f ", G(a, imaxLocal + 1, x, y));
                break;
            case BOTTOM:
                fprintf(fh, "%12.8f ", G(a, x, 0, y));
                break;
            case TOP:
                fprintf(fh, "%12.8f ", G(a, x, jmaxLocal + 1, y));
                break;
            case FRONT:
                fprintf(fh, "%12.8f ", G(a, x, y, 0));
                break;
            case BACK:
                fprintf(fh, "%12.8f ", G(a, x, y, kmaxLocal + 1));
                break;
            case NDIRS:
                printf("ERROR\n");
                break;
            }
        }
        fprintf(fh, "\n");
    }
    fclose(fh);
}

void testPrintHalo(Solver* s, double* a)
{
    int imaxLocal = s->comm->imaxLocal;
    int jmaxLocal = s->comm->jmaxLocal;
    int kmaxLocal = s->comm->kmaxLocal;

    printPlane(s, a, kmaxLocal + 2, imaxLocal + 2, BOTTOM);
    printPlane(s, a, kmaxLocal + 2, imaxLocal + 2, TOP);
    printPlane(s, a, kmaxLocal + 2, jmaxLocal + 2, LEFT);
    printPlane(s, a, kmaxLocal + 2, jmaxLocal + 2, RIGHT);
    printPlane(s, a, jmaxLocal + 2, imaxLocal + 2, FRONT);
    printPlane(s, a, jmaxLocal + 2, imaxLocal + 2, BACK);
}
