#ifndef MAIN_HEAD_H
#define MAIN_HEAD_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <stdarg.h>
#include <getopt.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <time.h>

// check if compiler understands OMP, if not, this file does probably not exist
#ifdef _OPENMP
    #include <omp.h>  
#endif

#define HDF5
#ifdef HDF5
    #include "hdf5.h"
#endif

#include "focal-struct.h"
#include "alloc-memory.h"
#include "constants.h"
#include "initialize_grid.h"
#include "background_profiles.h"
#include "advance_fields.h"
#include "boundary_module.h"
#include "save_data.h"
#include "free_allocated_memory.h"
#include "power_module.h"
#include "antenna_module.h"
#include "antenna_detector_module.h"

#endif