#ifndef ANTENNA_DETECTOR_MODULE_H
#define ANTENNA_DETECTOR_MODULE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "alloc-memory.h"
#include "focal-struct.h"
#include "macros-grid.h"
#include "constants.h"

void initialize_antDetect(  antennaDetector *antDetect, 
                            gridConfiguration *gridCfg, 
                            beamAntennaConfiguration *beamAnt);

void print_antDetect_info(antennaDetector *ant_Detect);

void init_antennaDetect(    antennaDetector *antDetect, 
                            gridConfiguration *gridCfg, 
                            beamAntennaConfiguration *beamAnt);

void control_antenna_detector( antennaDetector *antDetect, gridConfiguration *gridCfg, systemGrid *G, int t_int );

int detAnt_01_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int );

int detAnt_02_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int );

int detAnt_03_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int );

int detAnt_04_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int );


#endif