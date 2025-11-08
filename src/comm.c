/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "allocate.h"
#include "comm.h"

#if defined(_MPI)
// subroutines local to this module
static int sizeOfRank(int rank, int size, int N)
{
    return N / size + ((N % size > rank) ? 1 : 0);
}

static void setupCommunication(Comm *c, Direction direction, int layer)
{
    int imaxLocal  = c->imaxLocal;
    int jmaxLocal  = c->jmaxLocal;
    int kmaxLocal  = c->kmaxLocal;

    size_t dblsize = sizeof(double);
    int sizes[NDIMS];
    int subSizes[NDIMS];
    int starts[NDIMS];
    int offset  = 0;

    sizes[IDIM] = imaxLocal + 2;
    sizes[JDIM] = jmaxLocal + 2;
    sizes[KDIM] = kmaxLocal + 2;

    if (layer == HALO) {
        offset = 1;
    }

    switch (direction) {
    case LEFT:
        subSizes[IDIM] = 1;
        subSizes[JDIM] = jmaxLocal;
        subSizes[KDIM] = kmaxLocal;
        starts[IDIM]   = 1 - offset;
        starts[JDIM]   = 1;
        starts[KDIM]   = 1;
        break;
    case RIGHT:
        subSizes[IDIM] = 1;
        subSizes[JDIM] = jmaxLocal;
        subSizes[KDIM] = kmaxLocal;
        starts[IDIM]   = imaxLocal + offset;
        starts[JDIM]   = 1;
        starts[KDIM]   = 1;
        break;
    case BOTTOM:
        subSizes[IDIM] = imaxLocal;
        subSizes[JDIM] = 1;
        subSizes[KDIM] = kmaxLocal;
        starts[IDIM]   = 1;
        starts[JDIM]   = 1 - offset;
        starts[KDIM]   = 1;
        break;
    case TOP:
        subSizes[IDIM] = imaxLocal;
        subSizes[JDIM] = 1;
        subSizes[KDIM] = kmaxLocal;
        starts[IDIM]   = 1;
        starts[JDIM]   = jmaxLocal + offset;
        starts[KDIM]   = 1;
        break;
    case FRONT:
        subSizes[IDIM] = imaxLocal;
        subSizes[JDIM] = jmaxLocal;
        subSizes[KDIM] = 1;
        starts[IDIM]   = 1;
        starts[JDIM]   = 1;
        starts[KDIM]   = 1 - offset;
        break;
    case BACK:
        subSizes[IDIM] = imaxLocal;
        subSizes[JDIM] = jmaxLocal;
        subSizes[KDIM] = 1;
        starts[IDIM]   = 1;
        starts[JDIM]   = 1;
        starts[KDIM]   = kmaxLocal + offset;
        break;
    case NDIRS:
        printf("ERROR!\n");
        break;
    }

    if (layer == HALO) {
        MPI_Type_create_subarray(NDIMS,
            sizes,
            subSizes,
            starts,
            MPI_ORDER_C,
            MPI_DOUBLE,
            &c->rbufferTypes[direction]);
        MPI_Type_commit(&c->rbufferTypes[direction]);
    } else if (layer == BULK) {
        MPI_Type_create_subarray(NDIMS,
            sizes,
            subSizes,
            starts,
            MPI_ORDER_C,
            MPI_DOUBLE,
            &c->sbufferTypes[direction]);
        MPI_Type_commit(&c->sbufferTypes[direction]);
    }
}

static void assembleResult(Comm *c,
    double *src,
    double *dst,
    int imaxLocal[],
    int jmaxLocal[],
    int kmaxLocal[],
    int offset[],
    int kmax,
    int jmax,
    int imax)
{
    int numRequests = 1;

    if (c->rank == 0) {
        numRequests = c->size + 1;
    }

    MPI_Request requests[numRequests];

    /* all ranks send their interpolated bulk array */
    MPI_Isend(src,
        c->imaxLocal * c->jmaxLocal * c->kmaxLocal,
        MPI_DOUBLE,
        0,
        0,
        c->comm,
        &requests[0]);

    /* rank 0 assembles the subdomains */
    if (c->rank == 0) {
        for (int i = 0; i < c->size; i++) {
            MPI_Datatype domainType;
            int oldSizes[NDIMS] = { kmax, jmax, imax };
            int newSizes[NDIMS] = { kmaxLocal[i], jmaxLocal[i], imaxLocal[i] };
            int starts[NDIMS]   = { offset[i * NDIMS + KDIM],
                  offset[i * NDIMS + JDIM],
                  offset[i * NDIMS + IDIM] };
            MPI_Type_create_subarray(
                NDIMS, oldSizes, newSizes, starts, MPI_ORDER_C, MPI_DOUBLE, &domainType);
            MPI_Type_commit(&domainType);

            MPI_Irecv(dst, 1, domainType, i, 0, c->comm, &requests[i + 1]);
            MPI_Type_free(&domainType);
        }
    }

    MPI_Waitall(numRequests, requests, MPI_STATUSES_IGNORE);
}

// subroutines local to this module
static int sum(int *sizes, int init, int offset, int coord)
{
    int sum = 0;

    for (int i = init - offset; coord > 0; i -= offset, --coord) {
        sum += sizes[i];
    }

    return sum;
}
#endif // defined _MPI

// exported subroutines
void commReduction(double *v, int op)
{
#if defined(_MPI)
    if (op == MAX) {
        MPI_Allreduce(MPI_IN_PLACE, v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    } else if (op == SUM) {
        MPI_Allreduce(MPI_IN_PLACE, v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
}

int commIsBoundary(Comm *c, Direction direction)
{
#if defined(_MPI)
    switch (direction) {
    case LEFT:
        return c->coords[ICORD] == 0;
        break;
    case RIGHT:
        return c->coords[ICORD] == (c->dims[ICORD] - 1);
        break;
    case BOTTOM:
        return c->coords[JCORD] == 0;
        break;
    case TOP:
        return c->coords[JCORD] == (c->dims[JCORD] - 1);
        break;
    case FRONT:
        return c->coords[KCORD] == 0;
        break;
    case BACK:
        return c->coords[KCORD] == (c->dims[KCORD] - 1);
        break;
    case NDIRS:
        printf("ERROR!\n");
        break;
    }
#endif

    return 1;
}

void commExchange(Comm *c, double *grid)
{
#if defined(_MPI)
    int counts[6]      = { 1, 1, 1, 1, 1, 1 };
    MPI_Aint displs[6] = { 0, 0, 0, 0, 0, 0 };

    MPI_Neighbor_alltoallw(grid,
        counts,
        displs,
        c->sbufferTypes,
        grid,
        counts,
        displs,
        c->rbufferTypes,
        c->comm);
#endif
}

void commShift(Comm *c, double *f, double *g, double *h)
{
#if defined(_MPI)
    MPI_Request requests[6] = { MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL };

    /* shift G */
    /* receive ghost cells from bottom neighbor */
    MPI_Irecv(
        g, 1, c->rbufferTypes[BOTTOM], c->neighbours[BOTTOM], 0, c->comm, &requests[0]);

    /* send ghost cells to top neighbor */
    MPI_Isend(g, 1, c->sbufferTypes[TOP], c->neighbours[TOP], 0, c->comm, &requests[1]);

    /* shift F */
    /* receive ghost cells from left neighbor */
    MPI_Irecv(f, 1, c->rbufferTypes[LEFT], c->neighbours[LEFT], 1, c->comm, &requests[2]);

    /* send ghost cells to right neighbor */
    MPI_Isend(
        f, 1, c->sbufferTypes[RIGHT], c->neighbours[RIGHT], 1, c->comm, &requests[3]);

    /* shift H */
    /* receive ghost cells from front neighbor */
    MPI_Irecv(
        h, 1, c->rbufferTypes[FRONT], c->neighbours[FRONT], 2, c->comm, &requests[4]);

    /* send ghost cells to back neighbor */
    MPI_Isend(h, 1, c->sbufferTypes[BACK], c->neighbours[BACK], 2, c->comm, &requests[5]);

    MPI_Waitall(6, requests, MPI_STATUSES_IGNORE);
#endif
}

void commGetOffsets(Comm *c, int offsets[], int kmax, int jmax, int imax)
{
#if defined(_MPI)
    int sum = 0;

    for (int i = 0; i < c->coords[ICORD]; i++) {
        sum += sizeOfRank(i, c->dims[ICORD], imax);
    }
    offsets[IDIM] = sum;
    sum           = 0;

    for (int i = 0; i < c->coords[JCORD]; i++) {
        sum += sizeOfRank(i, c->dims[JCORD], jmax);
    }
    offsets[JDIM] = sum;
    sum           = 0;

    for (int i = 0; i < c->coords[KCORD]; i++) {
        sum += sizeOfRank(i, c->dims[KCORD], kmax);
    }
    offsets[KDIM] = sum;
#endif
}

#define G(v, i, j, k)                                                                    \
    v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

void commCollectResult(Comm *c,
    double *ug,
    double *vg,
    double *wg,
    double *pg,
    double *u,
    double *v,
    double *w,
    double *p,
    int kmax,
    int jmax,
    int imax)
{
    int imaxLocal = c->imaxLocal;
    int jmaxLocal = c->jmaxLocal;
    int kmaxLocal = c->kmaxLocal;
#if defined(_MPI)
    int offset[c->size * NDIMS];
    int imaxLocalAll[c->size];
    int jmaxLocalAll[c->size];
    int kmaxLocalAll[c->size];

    MPI_Gather(&imaxLocal, 1, MPI_INT, imaxLocalAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&jmaxLocal, 1, MPI_INT, jmaxLocalAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&kmaxLocal, 1, MPI_INT, kmaxLocalAll, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (c->rank == 0) {
        for (int i = 0; i < c->size; i++) {
            int coords[NCORDS];
            MPI_Cart_coords(c->comm, i, NDIMS, coords);
            offset[i * NDIMS + IDIM] =
                sum(imaxLocalAll, i, c->dims[IDIM] * c->dims[JDIM], coords[ICORD]);
            offset[i * NDIMS + JDIM] = sum(jmaxLocalAll, i, c->dims[IDIM], coords[JCORD]);
            offset[i * NDIMS + KDIM] = sum(kmaxLocalAll, i, 1, coords[KCORD]);

            printf("Rank: %d, Coords(k,j,i): %d %d %d, Size(k,j,i): %d %d %d, "
                   "Offset(k,j,i): %d %d %d\n",
                i,
                coords[KCORD],
                coords[JCORD],
                coords[ICORD],
                kmaxLocalAll[i],
                jmaxLocalAll[i],
                imaxLocalAll[i],
                offset[i * NDIMS + KDIM],
                offset[i * NDIMS + JDIM],
                offset[i * NDIMS + IDIM]);
        }
    }

    size_t bytesize = imaxLocal * jmaxLocal * kmaxLocal * sizeof(double);
    double *tmp     = allocate(64, bytesize);
    int idx         = 0;

    /* collect P */
    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                tmp[idx++] = G(p, i, j, k);
            }
        }
    }

    assembleResult(
        c, tmp, pg, imaxLocalAll, jmaxLocalAll, kmaxLocalAll, offset, kmax, jmax, imax);

    /* collect U */
    idx = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                tmp[idx++] = (G(u, i, j, k) + G(u, i - 1, j, k)) / 2.0;
            }
        }
    }

    assembleResult(
        c, tmp, ug, imaxLocalAll, jmaxLocalAll, kmaxLocalAll, offset, kmax, jmax, imax);

    /* collect V */
    idx = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                tmp[idx++] = (G(v, i, j, k) + G(v, i, j - 1, k)) / 2.0;
            }
        }
    }

    assembleResult(
        c, tmp, vg, imaxLocalAll, jmaxLocalAll, kmaxLocalAll, offset, kmax, jmax, imax);

    /* collect W */
    idx = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                tmp[idx++] = (G(w, i, j, k) + G(w, i, j, k - 1)) / 2.0;
            }
        }
    }

    assembleResult(
        c, tmp, wg, imaxLocalAll, jmaxLocalAll, kmaxLocalAll, offset, kmax, jmax, imax);

    free(tmp);
#else
    int idx = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                pg[idx++] = G(p, i, j, k);
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                ug[idx++] = (G(u, i, j, k) + G(u, i - 1, j, k)) / 2.0;
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                vg[idx++] = (G(v, i, j, k) + G(v, i, j - 1, k)) / 2.0;
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                wg[idx++] = (G(w, i, j, k) + G(w, i, j, k - 1)) / 2.0;
            }
        }
    }
#endif
}

void commPrintConfig(Comm *c)
{
#if defined(_MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    if (commIsMaster(c)) {
        printf("Communication setup:\n");
    }

    for (int i = 0; i < c->size; i++) {
        if (i == c->rank) {
            printf("\tRank %d of %d\n", c->rank, c->size);
            printf("\tNeighbours (front, back, bottom, top, left, right): %d, "
                   "%d, %d, "
                   "%d, %d, %d\n",
                c->neighbours[FRONT],
                c->neighbours[BACK],
                c->neighbours[BOTTOM],
                c->neighbours[TOP],
                c->neighbours[LEFT],
                c->neighbours[RIGHT]);
            printf("\tCoordinates (k,j,i) %d %d %d\n",
                c->coords[KCORD],
                c->coords[JCORD],
                c->coords[ICORD]);
            printf("\tLocal domain size (k,j,i) %dx%dx%d\n",
                c->kmaxLocal,
                c->jmaxLocal,
                c->imaxLocal);
            fflush(stdout);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void commInit(Comm *c, int argc, char **argv)
{
#if defined(_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(c->rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(c->size));
#else
    c->rank = 0;
    c->size = 1;
#endif
}

void commPartition(Comm *c, int kmax, int jmax, int imax)
{
#if defined(_MPI)
    int dims[NDIMS]    = { 0, 0, 0 };
    int periods[NDIMS] = { 0, 0, 0 };
    MPI_Dims_create(c->size, NDIMS, dims);
    MPI_Cart_create(MPI_COMM_WORLD, NCORDS, dims, periods, 0, &c->comm);
    MPI_Cart_shift(c->comm, ICORD, 1, &c->neighbours[LEFT], &c->neighbours[RIGHT]);
    MPI_Cart_shift(c->comm, JCORD, 1, &c->neighbours[BOTTOM], &c->neighbours[TOP]);
    MPI_Cart_shift(c->comm, KCORD, 1, &c->neighbours[FRONT], &c->neighbours[BACK]);
    MPI_Cart_get(c->comm, NCORDS, c->dims, periods, c->coords);

    c->imaxLocal = sizeOfRank(c->coords[KDIM], dims[ICORD], imax);
    c->jmaxLocal = sizeOfRank(c->coords[JDIM], dims[JCORD], jmax);
    c->kmaxLocal = sizeOfRank(c->coords[IDIM], dims[KCORD], kmax);

    // setup buffer types for communication
    setupCommunication(c, LEFT, BULK);
    setupCommunication(c, LEFT, HALO);
    setupCommunication(c, RIGHT, BULK);
    setupCommunication(c, RIGHT, HALO);
    setupCommunication(c, BOTTOM, BULK);
    setupCommunication(c, BOTTOM, HALO);
    setupCommunication(c, TOP, BULK);
    setupCommunication(c, TOP, HALO);
    setupCommunication(c, FRONT, BULK);
    setupCommunication(c, FRONT, HALO);
    setupCommunication(c, BACK, BULK);
    setupCommunication(c, BACK, HALO);
#else
    c->imaxLocal = imax;
    c->jmaxLocal = jmax;
    c->kmaxLocal = kmax;
#endif
}

void commFinalize(Comm *c)
{
#if defined(_MPI)
    for (int i = 0; i < NDIRS; i++) {
        MPI_Type_free(&c->sbufferTypes[i]);
        MPI_Type_free(&c->rbufferTypes[i]);
    }

    MPI_Finalize();
#endif
}

void commUpdateDatatypes(
    Comm *oldcomm, Comm *newcomm, int imaxLocal, int jmaxLocal, int kmaxLocal)
{
#if defined _MPI

    int result = MPI_Comm_dup(oldcomm->comm, &newcomm->comm);

    if (result == MPI_ERR_COMM) {
        printf("\nNull communicator. Duplication failed !!\n");
    }

    newcomm->rank      = oldcomm->rank;
    newcomm->size      = oldcomm->size;

    newcomm->imaxLocal = imaxLocal / 2;
    newcomm->jmaxLocal = jmaxLocal / 2;
    newcomm->kmaxLocal = kmaxLocal / 2;

    setupCommunication(newcomm, LEFT, BULK);
    setupCommunication(newcomm, LEFT, HALO);
    setupCommunication(newcomm, RIGHT, BULK);
    setupCommunication(newcomm, RIGHT, HALO);
    setupCommunication(newcomm, BOTTOM, BULK);
    setupCommunication(newcomm, BOTTOM, HALO);
    setupCommunication(newcomm, TOP, BULK);
    setupCommunication(newcomm, TOP, HALO);
    setupCommunication(newcomm, FRONT, BULK);
    setupCommunication(newcomm, FRONT, HALO);
    setupCommunication(newcomm, BACK, BULK);
    setupCommunication(newcomm, BACK, HALO);
#else
    newcomm->imaxLocal = imaxLocal;
    newcomm->jmaxLocal = jmaxLocal;
    newcomm->kmaxLocal = kmaxLocal;
#endif
}

void commFreeCommunicator(Comm *comm)
{
#ifdef _MPI
    MPI_Comm_free(&comm->comm);
#endif
}
