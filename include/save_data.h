#ifndef SAVE_DATA_H
#define SAVE_DATA_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "grid_io.h"
#include "auxiliar_module.h"

void simulation_folder(const char *path);
void system_evolution_folder(const char *path, const char *folder_name);
void data_folder(const char *path, const char *folder_name);
void copyJSON(const char *path, const char *folder_name);
void create_folder(saveData *saveDCfg);

void save_data_Grid( gridConfiguration *gridCfg, saveData *saveDCfg, systemGrid *G, int t_int );

void save_data_Configuration(   gridConfiguration *gridCfg, 
                                saveData *saveDCfg, 
                                systemGrid *G, 
                                beamAntennaConfiguration *beamAnt );

#endif