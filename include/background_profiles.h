#ifndef BACKGROUND_PROFILES_H
#define BACKGROUND_PROFILES_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "hdf5.h"
#include "grid_io.h"
#include "focal-struct.h"
#include "macros-grid.h"


int make_density_profile( gridConfiguration *gridCfg,
                          systemGrid *G,
                          double cntrl_para);
                          
int make_B0_profile( gridConfiguration *gridCfg, 
                     systemGrid *G,
                     double cntrl_para );

void control_background_profiles(gridConfiguration *gridCfg,
                             systemGrid *G);

#endif