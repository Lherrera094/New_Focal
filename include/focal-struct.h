#ifndef FOCAL_STRUCT_H
#define FOCAL_STRUCT_H

// define structures
typedef struct gridConfiguration {      /*Saves the main physical parameters of the system*/
    int
        Nx, Ny, Nz,   
        Nz_ref,
        d_absorb,
        t_end,
        ne_profile, B0_profile,
        boundary;
    double
        period,
        dx,dt,
        ne_0, B0_value;
} gridConfiguration;

typedef struct systemGrid{              /*Stores the wave components, plasma current, magnetic field and plasma density*/
    double *EB_WAVE, *EB_WAVE_ref,
           *J_B0,    *n_e;
} systemGrid;

typedef struct saveData{                /*Variables related to simulation data saving*/
    double *data2save,
           *timetraces;
    const char 	*foldername, *projectPath,
    		    *file_hdf5, *file_trace,
                *file_config;
    int     t_save;
} saveData;

typedef struct boundaryGrid{            /*Store grid value for the boundary variables*/

    /*ABC boundary*/
    double eco,
    
    /*Mur boundary*/
    *E_Xdir_OLD, *E_Ydir_OLD, *E_Zdir_OLD,
    *E_Xdir_OLD_ref, *E_Ydir_OLD_ref, *E_Zdir_OLD_ref,

    /*UPML boundary*/
    *DH_WAVE, *DH_WAVE_ref,
    *F1x, *F1y, *F1z,
    *F2x, *F2y, *F2z,
    *Cx, *Cy, *Cz,
    *F1zr, *F2zr,
    *Czr;
} boundaryGrid;

typedef struct beamAntennaConfiguration {   /*Antenna configuration variables*/
    int
        T_wave,
        exc_signal,
        ant_x, ant_y, ant_z,
        rampUpMethod;
    double
        omega_t,
        *antField_xy, *antPhaseTerms,
        antAngle_zy, antAngle_zx,
        ant_w0x, ant_w0y,
        z2waist;
} beamAntennaConfiguration;

typedef struct antennaDetector{
    int 
        antDetect_1D, 
        detAnt_01_zpos, detAnt_02_zpos,
        detAnt_03_zpos, detAnt_04_zpos,
        detAnt_01_ypos;
    double
        *detAnt_01_fields,
        *detAnt_02_fields,
        *detAnt_03_fields,
        *detAnt_04_fields;
} antennaDetector;

typedef struct powerCalcValues{
    int 
        pwr_dect;

    double
        power_abs_x1, power_abs_x2,
        power_abs_y1, power_abs_y2,
        power_abs_z1, power_abs_z2,
        power_abs_ref,

        poynt_x1, poynt_x2,
        poynt_y1, poynt_y2,
        poynt_z1, poynt_z2,
        poynt_z1_ref;    
} powerCalcValues;

typedef struct codeDiagnostics{

    double 
        *E, *R,
        *A_ref, *A_inc,
        *PML_values;
} codeDiagnostics;

enum FieldID {
    FIELD_01,
    FIELD_02,
    FIELD_03,
    FIELD_04
};

#endif
