// header guard first to prevent multiple declarations
#ifndef PHYSICS_MODULE_H
#define PHYSICS_MODULE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "constants.h"


int advance_J( gridConfiguration *gridCfg, 
               systemGrid *G );

int advance_B( gridConfiguration *gridCfg, 
               systemGrid *G );

int advance_B_ref( gridConfiguration *gridCfg, 
                   systemGrid *G);

int advance_E( gridConfiguration *gridCfg, 
               systemGrid *G );

int advance_E_ref( gridConfiguration *gridCfg, 
                   systemGrid *G ); 

void advance_fields( gridConfiguration *gridCfg, 
                    systemGrid *G );

int set2zero_1D( size_t N_x, double arr_1D[N_x] );
int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] );

#endif  // FOCAL_H

