#include "main_head.h"

int main( int argc, char *argv[] ){

    clock_t start, end;
    double cpu_time_used;
    start = clock();
    
#ifdef _OPENMP
    int n_threads;                          // number of threads that will be used (OpenMP)
#endif

    /*Call structs*/
    gridConfiguration           gridCfg;
    powerCalcValues             powerValStr; 
    systemGrid                  *G;
    saveData                    *saveDCfg;   
    boundaryGrid                *boundaryG;
    beamAntennaConfiguration    *beamAnt;
    antennaDetector             *antDetect;

    /*Alloc structs in memory*/
    ALLOC_1D( G, 1, systemGrid);
    ALLOC_1D( saveDCfg, 1, saveData);
    ALLOC_1D( boundaryG, 1, boundaryGrid);
    ALLOC_1D( beamAnt, 1, beamAntennaConfiguration);
    ALLOC_1D( antDetect, 1, antennaDetector );
    
    /*Initialize system values*/
    control_gridInit( &gridCfg, G, saveDCfg, boundaryG, beamAnt, antDetect);/*Initialize grid system matrices*/
    init_boundary(    &gridCfg, boundaryG);                                 /*Initialize Boundary conditions*/
    create_folder(    saveDCfg);                                            /*Call function to create folders*/
    init_antenna(     &gridCfg, beamAnt );                                  /*Initialize injection antenna*/
    initialize_power( &gridCfg, &powerValStr );                             /*Init power values*/
    control_background_profiles( &gridCfg, G);                              /*Define background ne and B0 profiles*/

    print_systemConfiguration( &gridCfg, beamAnt);                          /*Print all system configurations*/
    init_antennaDetect( antDetect, &gridCfg, beamAnt);

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
        
        wave_injection( &gridCfg, beamAnt, G, t_int );                          /*antenna_module.c*/
        advance_boundary(&gridCfg, G, boundaryG);                               /*boundary_module.c*/
        advance_fields(&gridCfg, G, boundaryG);                                 /*advance_fields.c*/

        control_antenna_detector( antDetect, &gridCfg, G, t_int );              /*antenna_detector_module.c*/
        control_power( &gridCfg, G, &powerValStr, saveDCfg, beamAnt, t_int );   /*power_module.c*/
        save_data_Grid( &gridCfg, saveDCfg, G, t_int );                         /*save_data.c*/
    }

    /*Save simulation data*/
    control_writeHDF5( &gridCfg, saveDCfg, G, beamAnt,antDetect );

    /*Free allocated memmory*/
    free_allocated_memory( &gridCfg, saveDCfg, G , boundaryG, beamAnt);

    end = clock();
    cpu_time_used = ( (double)(end-start) )/CLOCKS_PER_SEC;
	printf("Program Running Time = %.2e s.\n", cpu_time_used);
    
    return EXIT_SUCCESS;
}
