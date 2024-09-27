#ifndef POWER_MODULE_H
#define POWER_MODULE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "focal-struct.h"
#include "macros-grid.h"

void initialize_power( gridConfiguration *gridCfg, powerCalcValues *powerValStr );

void control_power(     gridConfiguration *gridCfg, 
                        systemGrid *G, 
                        powerCalcValues *powerValStr,
                        saveData *saveDCfg, 
                        beamAntennaConfiguration *beamAnt,
                        int t_int );
void compute_power( gridConfiguration *gridCfg, systemGrid *G, powerCalcValues *powerValStr, int t_int );
double calc_poynt_4(    gridConfiguration *gridCfg,
                        systemGrid *G,
                        powerCalcValues *powerValStr,
                        char absorber[] );

void compute_timetraces(    gridConfiguration *gridCfg, 
                            saveData *saveDCfg, 
                            beamAntennaConfiguration *beamAnt,
                            powerCalcValues *powerValStr,
                            int t_int );

#endif