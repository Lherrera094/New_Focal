#ifndef AUXILIAR_MODULE_H
#define AUXILIAR_MODULE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "focal-struct.h"
#include "macros-grid.h"

int set2zero_1D( size_t N_x, double arr_1D[N_x] );
int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] );
int setZero2save( gridConfiguration *gridCfg, saveData *saveDCfg);

#endif