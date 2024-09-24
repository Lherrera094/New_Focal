#ifndef POWER_MODULE_H
#define POWER_MODULE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "focal-struct.h"
#include "macros-grid.h"

void compute_timetraces( gridConfiguration *gridCfg, saveData *saveDCfg, int t_int );

#endif