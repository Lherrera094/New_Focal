#ifndef GRID_IO_H
#define GRID_IO_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "hdf5.h"
#include "focal-struct.h"
#include "macros-grid.h"

/*Read external density profile*/
int read_ProfileHDF( char filename[], char dataset[], double *array_3D );

/*Write Timetraces*/
int writeTimetraces2ascii( gridConfiguration *gridCfg, saveData *saveDCfg, char fullDirectory[] );

/*Write hdf5 files*/
int writeConfig2HDF( gridConfiguration *gridCfg, saveData *saveDCfg, char filename[], beamAntennaConfiguration *beamAnt );
int writeMyHDF_v4( gridConfiguration *gridCfg, char filename[], char dataset[], double *array_3D );

#endif