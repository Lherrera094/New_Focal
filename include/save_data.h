#ifndef SAVE_DATA_H
#define SAVE_DATA_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "grid_io.h"
#include "auxiliar_module.h"

void create_folder(saveData *saveDCfg);
void simulation_folder(const char *path);
void system_evolution_folder(const char *path, const char *folder_name);
void data_folder(const char *path, const char *folder_name);
void copyJSON(const char *path, const char *folder_name);

void control_writeHDF5( gridConfiguration *gridCfg, 
                        saveData *saveDCfg, 
                        systemGrid *G, 
                        beamAntennaConfiguration *beamAnt,
                        antennaDetector *antDetect );

void save_data_Grid( gridConfiguration *gridCfg, saveData *saveDCfg, systemGrid *G, int t_int );

void save_data_Configuration(   gridConfiguration *gridCfg, 
                                saveData *saveDCfg, 
                                systemGrid *G, 
                                beamAntennaConfiguration *beamAnt,
                                char filename[] );

void save_antennaDetect(    gridConfiguration *gridCfg,
                            antennaDetector *antDetect,
                            char filename_hdf5[]);

void writeFile( gridConfiguration *gridCfg, codeDiagnostics *diagnostict, saveData *saveDCfg );

void writeUPMLdata( gridConfiguration *gridCfg, codeDiagnostics *diagnostic, 
                    saveData *saveDCfg, systemGrid *G, int t_int );

#endif