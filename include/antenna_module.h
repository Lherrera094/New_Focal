#ifndef ANTENNA_MODULE_H
#define ANTENNA_MODULE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "constants.h"

/*Initialize antenna*/
void init_antenna( gridConfiguration *gridCfg, beamAntennaConfiguration *beamAnt );
int make_antenna_profile( gridConfiguration *gridCfg, beamAntennaConfiguration *beamAnt );

/*Inject wave*/
void wave_injection(    gridConfiguration *gridCfg, 
                        beamAntennaConfiguration *beamAnt, 
                        systemGrid *G, int t_int );
                        
int add_source( gridConfiguration *gridCfg, beamAntennaConfiguration *beamAnt, 
                systemGrid *G,
                double Y, int t_int );

int add_source_ref( gridConfiguration *gridCfg, beamAntennaConfiguration *beamAnt, 
                    systemGrid *G,
                    double Y, int t_int );

/*Auxiliar functions*/
double antenna_field_rampup( int rampUpMth, double wave_period, int t_int );

double antenna_calcHansenExEy_O( double theta_rad, double Y );

#endif