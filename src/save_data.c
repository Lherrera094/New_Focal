#include "save_data.h"

/*Call functions to create path folder and simulation folder*/
void simulation_folder(const char *path){
    
    struct stat st = {0};

    /*Checks if directory exists*/
    if( stat(path, &st) == -1){
        //Directory does not exists. Create it.
        if(mkdir(path, 0700) == 0){
            printf("Main project folder created successfully.\n");
        }else{
            printf("Error creating directory: %s\n", path);
            return;
        }
    }
}

void system_evolution_folder(const char *path, const char *folder_name){
    
    char fullPath[2024];

    // Create the full directory path and check for buffer overflow
    if (snprintf(fullPath, sizeof(fullPath), "%s/%s/System_evolution", path, folder_name) >= sizeof(fullPath)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return;
    }

    struct stat st = {0};

    /*Checks if directory exists.*/
    if( stat(fullPath, &st) == -1){
        //Directory does not exists. Create it.
        if( mkdir(fullPath, 0700) == 0){
            //printf("%s folder created successfully. \n", folder_name);
        }else{
            printf("Error creating directory: %s\n", folder_name);
            return;
        }
    }

}

void data_folder(const char *path, const char *folder_name){
    
    char fullPath[1024];

    // Create the full directory path and check for buffer overflow
    if (snprintf(fullPath, sizeof(fullPath), "%s/%s", path, folder_name) >= sizeof(fullPath)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return;
    }

    struct stat st = {0};

    /*Checks if directory exists.*/
    if( stat(fullPath, &st) == -1){
        //Directory does not exists. Create it.
        if( mkdir(fullPath, 0700) == 0){
            printf("%s folder created successfully. \n", folder_name);
        }else{
            printf("Error creating directory: %s\n", folder_name);
            return;
        }
    }else{
        printf("%s already exists.\n", folder_name);
    }
}

void copyJSON(const char *path, const char *folder_name){
    
    char destination[1024];

    //Read the source file
    FILE *srcFile = fopen("../input_FOCAL.json", "rb");
    if (srcFile == NULL) {
        perror("Error opening source file");
        return;
    }

    // Create the full directory path and check for buffer overflow
    if (snprintf(destination, sizeof(destination), "%s/%s/input_FOCAL.json", path, folder_name) >= sizeof(destination)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return;
    }

    //open destination file
    FILE *destFile = fopen(destination, "wb");
    if (destFile == NULL) {
        perror("Error opening destination file");
        fclose(srcFile);
        return;
    }

    char buffer[1024];
    size_t bytesRead;
    while ((bytesRead = fread(buffer, 1, sizeof(buffer), srcFile)) > 0) {
        fwrite(buffer, 1, bytesRead, destFile);
    }

    fclose(srcFile);
    fclose(destFile);

    printf("JSON file saved.\n");
}

void create_folder(saveData *saveDCfg){

    simulation_folder( projectPath );
    data_folder( projectPath, foldername );
    system_evolution_folder(projectPath, foldername );
    copyJSON( projectPath, foldername );
    
}

/*Functions to save data in folders*/
void control_writeHDF5( gridConfiguration *gridCfg, 
                        saveData *saveDCfg, 
                        systemGrid *G, 
                        beamAntennaConfiguration *beamAnt,
                        antennaDetector *antDetect ){

    /*Char values as directions to the correct folder*/
    char fullDir[1000], filename_config[2024];
    if (snprintf(fullDir, sizeof(fullDir),"%s/%s", projectPath, foldername) >= sizeof(fullDir)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
    }
    //Append the name of the files
    snprintf( filename_config, sizeof(filename_config), "%s/%s", fullDir, file_config);

    writeTimetraces2ascii( gridCfg, saveDCfg, fullDir );
    save_data_Configuration( gridCfg, saveDCfg, G, beamAnt, filename_config );
    save_antennaDetect( gridCfg, antDetect, filename_config );

}

void save_data_Configuration(   gridConfiguration *gridCfg, 
                                saveData *saveDCfg, 
                                systemGrid *G, 
                                beamAntennaConfiguration *beamAnt,
                                char filename[] ){
    
    size_t ii, jj, kk;

    /*Call functions to save simulation data*/ 
    writeConfig2HDF( gridCfg, saveDCfg, filename, beamAnt );

    /*Background Density*/
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg, filename, "n_e", G->n_e) );

    // background magnetic field
    // B0x: even-odd-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii < Nx ; ii+=2) {
        for (jj=0 ; jj < Ny ; jj+=2) {
            for (kk=0 ; kk < Nz ; kk+=2) {
                data2save(ii/2,jj/2,kk/2) = J_B0(ii  ,jj+1,kk+1);
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg, filename, "B0x", saveDCfg->data2save) );
    setZero2save( gridCfg, saveDCfg);

    // B0y: odd-even-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii < Nx ; ii+=2) {
        for (jj=0 ; jj < Ny ; jj+=2) {
            for (kk=0 ; kk < Nz ; kk+=2) {
                data2save(ii/2,jj/2,kk/2) = J_B0(ii+1,jj  ,kk+1);
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg, filename, "B0y", saveDCfg->data2save) );
    setZero2save( gridCfg, saveDCfg);

    // B0z: odd-odd-even
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii < Nx ; ii+=2) {
        for (jj=0 ; jj < Ny ; jj+=2) {
            for (kk=0 ; kk < Nz ; kk+=2) {
                data2save(ii/2,jj/2,kk/2) = J_B0(ii+1,jj+1,kk  );
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg, filename, "B0z", saveDCfg->data2save) ) ;
    setZero2save( gridCfg, saveDCfg);
}

void save_antennaDetect(    gridConfiguration *gridCfg,
                            antennaDetector *antDetect,
                            char filename_hdf5[]){

    if ( detAnt_01_z < ( Nz - d_absorb )) {
        detAnt1D_write2hdf5( Nx, antDetect, filename_hdf5, "/detAnt_01" , 
                             detAnt_01_y, detAnt_01_z,
                             FIELD_01 );
    }
    if ( detAnt_02_z < ( Nz - d_absorb )) {
        detAnt1D_write2hdf5( Nx, antDetect, filename_hdf5, "/detAnt_02" , 
                             detAnt_01_y, detAnt_02_z,
                             FIELD_02 );
    }
    if ( detAnt_03_z < ( Nz - d_absorb )) {
        detAnt1D_write2hdf5( Nx, antDetect, filename_hdf5, "/detAnt_03" , 
                             detAnt_01_y, detAnt_03_z,
                             FIELD_03 );
    }
    if ( detAnt_04_z < ( Nz - d_absorb )) {
        detAnt1D_write2hdf5( Nx, antDetect, filename_hdf5, "/detAnt_04" , 
                             detAnt_01_y, detAnt_04_z,
                             FIELD_04 );
    }

}

void save_data_Grid(    gridConfiguration *gridCfg, 
                        saveData *saveDCfg, 
                        systemGrid *G , int t_int){

    if( t_int % (t_save * (int)period) == 0 && t_int != 0){

        size_t ii, jj, kk;
        /*Char values as directions to the correct folder*/
        char fullDir[1000], filename_physics[2024];
        if (snprintf(fullDir, sizeof(fullDir),"%s/%s", projectPath, foldername) >= sizeof(fullDir)) {
            fprintf(stderr, "Error: Directory path is too long.\n");
        }
        //Append the name of the files
        snprintf( filename_physics, sizeof(filename_physics), "%s/System_evolution/%s_time=%d", fullDir, file_hdf5, t_int / (int)period );

        // save into hdf5
        // abs(E)
        // prepare array for that
    #pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
        for (ii=0 ; ii < Nx ; ii+=2) {
            for (jj=0 ; jj < Ny ; jj+=2) {
                for (kk=0 ; kk < Nz ; kk+=2) {
                    data2save(ii/2,jj/2,kk/2) = 
                        sqrt(  pow( EB_WAVE(ii+1,jj  ,kk  ), 2) 
                              +pow( EB_WAVE(ii  ,jj+1,kk  ), 2) 
                              +pow( EB_WAVE(ii  ,jj  ,kk+1), 2) );
                }
            }
        }
        printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg, filename_physics, "E_abs", saveDCfg->data2save) ) ;
        setZero2save( gridCfg, saveDCfg);

        // abs(B)
        // prepare array for that
    /*#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
        for (ii=0 ; ii < Nx ; ii+=2) {
            for (jj=0 ; jj < Ny ; jj+=2) {
                for (kk=0 ; kk < Nz ; kk+=2) {
                    data2save(ii/2,jj/2,kk/2) = 
                        sqrt(  pow( EB_WAVE(ii  ,jj+1,kk+1), 2) 
                              +pow( EB_WAVE(ii+1,jj  ,kk+1), 2) 
                              +pow( EB_WAVE(ii+1,jj+1,kk  ), 2) );
                }
            }
        }
        printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg, filename_physics, "B_abs", saveDCfg->data2save) ) ;
        */

    }//end if

}

