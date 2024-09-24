#include "free_allocated_memory.h"

// Function to free the allocated 3D array memory

void free_systemGrid(systemGrid *G){
    free(G->EB_WAVE);
    free(G->EB_WAVE_ref);
    free(G->J_B0);
    free(G->n_e);
}

void free_boundaryGrid( gridConfiguration *gridCfg, boundaryGrid *boundaryG ){

    if (boundary_sel == 2){

        free(boundaryG->E_Xdir_OLD);
        free(boundaryG->E_Ydir_OLD);
        free(boundaryG->E_Zdir_OLD);
        free(boundaryG->E_Xdir_OLD_ref);
        free(boundaryG->E_Ydir_OLD_ref);
        free(boundaryG->E_Zdir_OLD_ref);

    }
    else if(boundary_sel == 3){
        
    }

}

void free_saveDCfg( saveData *saveDCfg ){

    free( saveDCfg->data2save );
    free( saveDCfg->timetraces );
    free( (void *)projectPath );
    free( (void *)foldername );
    free( (void *)file_hdf5 );
    free( (void *)file_trace );
    free( (void *)file_config );

}

void free_beamAnt( beamAntennaConfiguration *beamAnt ){

    free( beamAnt->antField_xy );
    free( beamAnt->antPhaseTerms );

}

void free_allocated_memory( gridConfiguration *gridCfg, 
                            saveData *saveDCfg, 
                            systemGrid *G,
                            boundaryGrid *boundaryG, 
                            beamAntennaConfiguration *beamAnt ){

    free_systemGrid( G );
    free_beamAnt( beamAnt );
    free_boundaryGrid( gridCfg, boundaryG );
    free(saveDCfg);
    free(beamAnt);
    free(G);

    printf("Allocated memory free. \n");
}