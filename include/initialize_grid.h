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
void gridConfInit(  gridConfiguration *gridCfg, 
                    saveData *saveDCfg, 
                    beamAntennaConfiguration *beamAnt,
                    antennaDetector *antDetect );

char *read_json();

void write_JSON_onStruct(   gridConfiguration *gridCfg, 
                            saveData *saveDCfg, 
                            beamAntennaConfiguration *beamAnt,
                            antennaDetector *antDetect );

void control_gridInit(  gridConfiguration *gridCfg, 
                        systemGrid *G, 
                        saveData *saveDCfg, 
                        boundaryGrid *boundaryG,
                        beamAntennaConfiguration *beamAnt,
                        antennaDetector *antDetect );

void print_systemConfiguration(gridConfiguration *gridCfg, beamAntennaConfiguration *beamAnt );

void allocate_data2save( gridConfiguration *gridCfg, saveData *saveDCfg);
void allocate_timetraces( gridConfiguration *gridCfg, saveData *saveDCfg );

int t_end_s(gridConfiguration *gridCfg);

#endif
