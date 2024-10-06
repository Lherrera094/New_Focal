#ifndef INITIALIZE_GRID_H
#define INITIALIZE_GRID_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cJSON.h"
#include "focal-struct.h"
#include "alloc-memory.h"
#include "macros-grid.h"

//functions in initialize grid
void control_gridInit(  gridConfiguration *gridCfg, 
                        systemGrid *G, 
                        saveData *saveDCfg, 
                        boundaryGrid *boundaryG,
                        beamAntennaConfiguration *beamAnt,
                        antennaDetector *antDetect );

void gridConfInit(  gridConfiguration *gridCfg, 
                    saveData *saveDCfg, 
                    beamAntennaConfiguration *beamAnt,
                    antennaDetector *antDetect );

char *read_json();

void write_JSON_onStruct(   gridConfiguration *gridCfg, 
                            saveData *saveDCfg, 
                            beamAntennaConfiguration *beamAnt,
                            antennaDetector *antDetect );

void allocateMemory_structs( gridConfiguration *gridCfg,
                             systemGrid *G, 
                             saveData *saveDCfg, 
                             beamAntennaConfiguration *beamAnt,
                             antennaDetector *antDetect );

void print_systemConfiguration(gridConfiguration *gridCfg, beamAntennaConfiguration *beamAnt );

int t_end_s(gridConfiguration *gridCfg);

#endif
