/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "comm.h"
#include "vtkWriter.h"

// reset fileview for output of string headers
static void resetFileview(VtkOptions *o)
{
    MPI_Offset disp;
    MPI_File_sync(o->fh);
    MPI_Barrier(o->comm.comm);
    MPI_File_get_size(o->fh, &disp);
    MPI_File_set_view(o->fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
}

static void writeVersion(VtkOptions *o)
{
    char header[50] = "# vtk DataFile Version 3.0\n";
    // always overwrite exiting files
    MPI_File_set_view(o->fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

    if (commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, header, (int)strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
    }
}

static void writeHeader(VtkOptions *o)
{
    char header[400];
    char *cursor = header;

    cursor += sprintf(cursor, "PAMPI cfd solver output\n");
    cursor += sprintf(cursor, "BINARY\n");

    cursor += sprintf(cursor, "DATASET STRUCTURED_POINTS\n");
    cursor += sprintf(
        cursor, "DIMENSIONS %d %d %d\n", o->grid->imax, o->grid->jmax, o->grid->kmax);
    cursor += sprintf(cursor,
        "ORIGIN %f %f %f\n",
        o->grid->dx * 0.5,
        o->grid->dy * 0.5,
        o->grid->dz * 0.5);
    cursor +=
        sprintf(cursor, "SPACING %f %f %f\n", o->grid->dx, o->grid->dy, o->grid->dz);
    cursor +=
        sprintf(cursor, "POINT_DATA %d\n", o->grid->imax * o->grid->jmax * o->grid->kmax);

    if (commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, header, (int)strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
    }
}

void vtkOpen(VtkOptions *o, char *problem)
{
    char filename[50];
    snprintf(filename, 50, "%s-p%d.vtk", problem, o->comm.size);
    MPI_File_open(
        o->comm.comm, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &o->fh);

    if (commIsMaster(&o->comm)) {
        printf("Writing VTK output for %s\n", problem);
    }

    writeVersion(o);
    writeHeader(o);
}

void vtkScalar(VtkOptions *o, char *name, double *s)
{
    resetFileview(o);
    if (commIsMaster(&o->comm))
        printf("Register scalar %s\n", name);

    char header[100];
    char *cursor = header;

    cursor += sprintf(cursor, "SCALARS %s double\n", name);

    if (commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, header, (int)strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    int offsets[NDIMS];
    commGetOffsets(&o->comm, offsets, o->grid->kmax, o->grid->jmax, o->grid->imax);

    // set global view in file
    MPI_Offset disp;
    MPI_Datatype fileViewType;
    MPI_File_sync(o->fh);
    MPI_Barrier(o->comm.comm);
    MPI_File_get_size(o->fh, &disp);

    MPI_Type_create_subarray(NDIMS,
        (int[NDIMS]) { o->grid->kmax, o->grid->jmax, o->grid->imax },
        (int[NDIMS]) { o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal },
        offsets,
        MPI_ORDER_C,
        MPI_DOUBLE,
        &fileViewType);
    MPI_Type_commit(&fileViewType);
    MPI_File_set_view(o->fh, disp, MPI_DOUBLE, fileViewType, "external32", MPI_INFO_NULL);

#ifdef VERBOSE
    printf("Rank: %d, Disp: %lld,  Size(k,j,i): %d %d %d, Offset(k,j,i): %d %d %d\n",
        o->comm.rank,
        disp,
        o->comm.kmaxLocal,
        o->comm.jmaxLocal,
        o->comm.imaxLocal,
        offsets[KDIM],
        offsets[JDIM],
        offsets[IDIM]);
#endif

    // create local bulk type
    MPI_Datatype bulkType;

    MPI_Type_create_subarray(NDIMS,
        (int[NDIMS]) { o->comm.kmaxLocal + 2,
            o->comm.jmaxLocal + 2,
            o->comm.imaxLocal + 2 }, // oldsizes
        (int[NDIMS]) {
            o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal }, // newsizes
        (int[NDIMS]) { 1, 1, 1 }, // offsets
        MPI_ORDER_C,
        MPI_DOUBLE,
        &bulkType);
    MPI_Type_commit(&bulkType);

    MPI_File_write(o->fh, s, 1, bulkType, MPI_STATUS_IGNORE);
    MPI_Type_free(&bulkType);
    MPI_Type_free(&fileViewType);

    // Binary segment must be terminated with newline character
    resetFileview(o);
    if (commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
    }
}

#define G(v, i, j, k)                                                                    \
    v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

void vtkVector(VtkOptions *o, char *name, VtkVector vec)
{
    int imaxLocal = o->comm.imaxLocal;
    int jmaxLocal = o->comm.jmaxLocal;
    int kmaxLocal = o->comm.kmaxLocal;

    if (commIsMaster(&o->comm))
        printf("Register vector %s\n", name);
    const size_t MAX_HEADER = 100;

    char *header            = (char *)malloc(MAX_HEADER);
    sprintf(header, "VECTORS %s double\n", name);

    resetFileview(o);
    if (commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, header, (int)strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    int offsets[NDIMS];
    commGetOffsets(&o->comm, offsets, o->grid->kmax, o->grid->jmax, o->grid->imax);

    // set global view in file
    MPI_Offset disp;
    MPI_Datatype fileViewType, vectorType;
    MPI_File_sync(o->fh);
    MPI_Barrier(o->comm.comm);
    MPI_File_get_size(o->fh, &disp);

    MPI_Type_contiguous(NDIMS, MPI_DOUBLE, &vectorType);
    MPI_Type_commit(&vectorType);

    MPI_Type_create_subarray(NDIMS,
        (int[NDIMS]) { o->grid->kmax, o->grid->jmax, o->grid->imax },
        (int[NDIMS]) { kmaxLocal, jmaxLocal, imaxLocal },
        offsets,
        MPI_ORDER_C,
        vectorType,
        &fileViewType);
    MPI_Type_commit(&fileViewType);
    MPI_File_set_view(o->fh, disp, MPI_DOUBLE, fileViewType, "external32", MPI_INFO_NULL);

    size_t cnt  = imaxLocal * jmaxLocal * kmaxLocal;
    double *tmp = allocate(64, cnt * NDIMS * sizeof(double));
    int idx     = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                tmp[idx++] = (G(vec.u, i, j, k) + G(vec.u, i - 1, j, k)) / 2.0;
                tmp[idx++] = (G(vec.v, i, j, k) + G(vec.v, i, j - 1, k)) / 2.0;
                tmp[idx++] = (G(vec.w, i, j, k) + G(vec.w, i, j, k - 1)) / 2.0;
            }
        }
    }

    if (commIsMaster(&o->comm))
        printf("Write %d vectors\n", (int)cnt);
    MPI_File_write(o->fh, tmp, cnt, vectorType, MPI_STATUS_IGNORE);
    MPI_Type_free(&fileViewType);
    MPI_Type_free(&vectorType);

    // Binary segment must be terminated with newline character
    resetFileview(o);
    if (commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
    }
}

void vtkClose(VtkOptions *o)
{
    MPI_File_close(&o->fh);
}
