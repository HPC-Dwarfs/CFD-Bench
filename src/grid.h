/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __GRID_H_
#define __GRID_H_

typedef struct {
    double dx, dy, dz;
    int imax, jmax, kmax;
    double xlength, ylength, zlength;
} Grid;

#endif // __GRID_H_
