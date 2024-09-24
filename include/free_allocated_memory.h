#ifndef FREE_ALLOCATED_MEMORY_H
#define FREE_ALLOCATED_MEMORY_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "focal-struct.h"
#include "macros-grid.h"

void free_systemGrid(systemGrid *G);
void free_data2save(gridConfiguration *gridCfg, saveData *saveDCfg);
void free_timetraces( gridConfiguration *gridCfg, saveData *saveDCfg );
void free_systemGrid(systemGrid *G);
void free_boundaryGrid( gridConfiguration *gridCfg, boundaryGrid *boundaryG );
void free_saveDCfg( saveData *saveDCfg );

void free_allocated_memory( gridConfiguration *gridCfg, 
                            saveData *saveDCfg, 
                            systemGrid *G,
                            boundaryGrid *boundaryG,
                            beamAntennaConfiguration *beamAnt );

#endif