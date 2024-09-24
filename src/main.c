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

int main( int argc, char *argv[] ){

    clock_t start, end;
    double cpu_time_used;
    start = clock();
    
#ifdef _OPENMP
    int n_threads;                          // number of threads that will be used (OpenMP)
#endif

    /*Call structs*/
    gridConfiguration           gridCfg;
    systemGrid                  *G;
    saveData                    *saveDCfg;   
    boundaryGrid                *boundaryG;
    beamAntennaConfiguration    *beamAnt;

    /*Alloc structs in memory*/
    ALLOC_1D( G, 1, systemGrid);
    ALLOC_1D( saveDCfg, 1, saveData);
    ALLOC_1D( boundaryG, 1, boundaryGrid);
    ALLOC_1D( beamAnt, 1, beamAntennaConfiguration);
    
    control_gridInit( &gridCfg, G, saveDCfg, boundaryG, beamAnt);            /*Initialize grid system matrices*/
    init_boundary(    &gridCfg, boundaryG);                                 /*Initialize Boundary conditions*/
    create_folder(    saveDCfg);                                            /*Call function to create folders*/
    init_antenna(     &gridCfg, beamAnt );                                  /*Initialize injection antenna*/
    control_background_profiles( &gridCfg, G);                              /*Define background ne and B0 profiles*/

    print_systemConfiguration( &gridCfg, beamAnt);                          /*Print all system configurations*/

#ifdef _OPENMP
#pragma omp parallel private(n_threads)
    {
    n_threads = omp_get_num_threads();
    #pragma omp single
    printf( "Number of threads that will be used (OpenMP) = %d\n", n_threads );
    }
#endif

    /*Initiate system evolution*/
    for (int t_int=0 ; t_int <= t_end_s(&gridCfg) ; ++t_int){
        
        wave_injection( &gridCfg, beamAnt, G, t_int );

        advance_boundary(&gridCfg, G, boundaryG);
        advance_fields(&gridCfg, G);
        
        compute_timetraces( &gridCfg, saveDCfg, t_int );
        save_data_Grid( &gridCfg, saveDCfg, G, t_int );
      
    }

    /*Save simulation data*/
    save_data_Configuration( &gridCfg, saveDCfg, G, beamAnt );

    /*Free allocated memmory*/
    free_allocated_memory( &gridCfg, saveDCfg, G , boundaryG, beamAnt);

    end = clock();
    cpu_time_used = ( (double)(end-start) )/CLOCKS_PER_SEC;
	printf("Program Running Time = %.2e s.\n", cpu_time_used);
    
    return EXIT_SUCCESS;

}
