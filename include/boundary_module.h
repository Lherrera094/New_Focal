#ifndef BOUNDARY_MODULE_H
#define BOUNDARY_MODULE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "focal-struct.h"
#include "alloc-memory.h"
#include "macros-grid.h"

#define ABSORBER_DAMPING(eco,damp) (1.-eco*damp*damp)

int apply_absorber( gridConfiguration *gridCfg, 
                    systemGrid *G, 
                    boundaryGrid *boundaryG );

int apply_absorber_ref( gridConfiguration *gridCfg, 
                        systemGrid *G, 
                        boundaryGrid *boundaryG );

int abc_Mur_saveOldE_xdir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG );

int abc_Mur_saveOldE_ydir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG );

int abc_Mur_saveOldE_zdir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG );

int abc_Mur_saveOldE_ref_xdir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG );

int abc_Mur_saveOldE_ref_ydir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG );

int abc_Mur_saveOldE_ref_zdir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG );


void init_boundary(gridConfiguration *gridCfg, boundaryGrid *boundaryG);
void advance_boundary(gridConfiguration *gridCfg, systemGrid *G, boundaryGrid *boundaryG);
int apply_numerical_viscosity( gridConfiguration *gridCfg, systemGrid *G );

#endif