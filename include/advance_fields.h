// header guard first to prevent multiple declarations
#ifndef PHYSICS_MODULE_H
#define PHYSICS_MODULE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "auxiliar_module.h"
#include "constants.h"
#include "UPML_module.h"

void advance_fields( gridConfiguration *gridCfg, 
                     systemGrid *G,
                     boundaryGrid *boundaryG );

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

int advance_J_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G );

int advance_B_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG );

int advance_Bref_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG );

int advance_E_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG );

int advance_Eref_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG );

#endif  // FOCAL_H

