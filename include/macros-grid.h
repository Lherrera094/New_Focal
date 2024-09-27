#ifndef MACROS_BOUNDARY_H
#define MACROS_BOUNDARY_H

/*Macros for grid configuration*/
#define NxG(gridCfg)            gridCfg->Nx
#define NyG(gridCfg)            gridCfg->Ny
#define NzG(gridCfg)            gridCfg->Nz
#define Nz_refG(gridCfg)        gridCfg->Nz_ref
#define d_absorbG(gridCfg)      gridCfg->d_absorb
#define t_endG(gridCfg)         gridCfg->t_end
#define periodG(gridCfg)        gridCfg->period
#define dxG(gridCfg)            gridCfg->dx
#define dtG(gridCfg)            gridCfg->dt
#define ne_profileG(gridCfg)    gridCfg->ne_profile
#define ne_0G(gridCfg)          gridCfg->ne_0
#define B0_profileG(gridCfg)    gridCfg->B0_profile
#define B0_valueG(gridCfg)      gridCfg->B0_value
#define boundaryG(gridCfg)      gridCfg->boundary

#define Nx                      NxG(gridCfg)             
#define Ny                      NyG(gridCfg)           
#define Nz                      NzG(gridCfg)  
#define Nz_ref                  Nz_refG(gridCfg)  
#define d_absorb                d_absorbG(gridCfg)
#define t_end                   t_endG(gridCfg)
#define period                  periodG(gridCfg)
#define dx                      dxG(gridCfg)
#define dt                      dtG(gridCfg)
#define ne_profile              ne_profileG(gridCfg)
#define ne_0                    ne_0G(gridCfg)
#define B0_profile              B0_profileG(gridCfg)
#define B0_value                B0_valueG(gridCfg)
#define boundary_sel            boundaryG(gridCfg)

/*Macros for Grid system*/
#define EB_WAVEg(G,i,j,k)             G->EB_WAVE[((i) * (Ny) + j) * (Nz) + k]
#define EB_WAVE_refg(G,i,j,k)         G->EB_WAVE_ref[((i) * (Ny) + j) * (Nz_ref) + k]
#define J_B0g(G,i,j,k)                G->J_B0[((i) * (Ny) + j) * (Nz) + k]
#define n_eg(G,i,j,k)                 G->n_e[((i) * (Ny/2) + j) * (Nz/2) + k]

#define EB_WAVE(i,j,k)                EB_WAVEg(G,i,j,k)
#define EB_WAVE_ref(i,j,k)            EB_WAVE_refg(G,i,j,k)
#define J_B0(i,j,k)                   J_B0g(G,i,j,k)
#define n_e(i,j,k)                    n_eg(G,i,j,k)

/*Macros for save data*/
#define data2saveSt(saveDCfg,i,j,k)     saveDCfg->data2save[((i) * (Ny/2) + j) * (Nz/2) + k]
#define timetracesSt(saveDCfg,i,j)      saveDCfg->timetraces[((i) * (8) ) + j]
#define projectPathSt(saveDCfg)         saveDCfg->projectPath
#define foldernameSt(saveDCfg)          saveDCfg->foldername
#define file_hdf5St(saveDCfg)           saveDCfg->file_hdf5
#define file_traceSt(saveDCfg)          saveDCfg->file_trace
#define file_ConfigSt(saveDCfg)         saveDCfg->file_config
#define t_saveSt(saveDCfg)              saveDCfg->t_save

#define data2save(i,j,k)                data2saveSt(saveDCfg,i,j,k)
#define timetraces(i,j)                 timetracesSt(saveDCfg,i,j)
#define projectPath                     projectPathSt(saveDCfg)
#define foldername                      foldernameSt(saveDCfg)
#define file_hdf5                       file_hdf5St(saveDCfg)
#define file_trace                      file_traceSt(saveDCfg)
#define file_config                     file_ConfigSt(saveDCfg)
#define t_save                          t_saveSt(saveDCfg)

/*Macros for antenna injection*/
#define antFieldBG(beamAnt,i,j)         beamAnt->antField_xy[( (i) * (Ny/2) ) + j ]
#define antPhaseBG(beamAnt,i,j)         beamAnt->antPhaseTerms[( (i) * (Ny/2) ) + j ]
#define T_waveBG(beamAnt)               beamAnt->T_wave
#define exc_signalBG(beamAnt)           beamAnt->exc_signal
#define ant_xBG(beamAnt)                beamAnt->ant_x
#define ant_yBG(beamAnt)                beamAnt->ant_y
#define ant_zBG(beamAnt)                beamAnt->ant_z
#define rampUpMBG(beamAnt)              beamAnt->rampUpMethod
#define omega_tBG(beamAnt)              beamAnt->omega_t
#define antAngle_zxBG(beamAnt)          beamAnt->antAngle_zx
#define antAngle_zyBG(beamAnt)          beamAnt->antAngle_zy
#define ant_w0xBG(beamAnt)              beamAnt->ant_w0x
#define ant_w0yBG(beamAnt)              beamAnt->ant_w0y
#define z2waistBG(beamAnt)              beamAnt->z2waist

#define antField_xy(i,j)                antFieldBG(beamAnt,i,j)         
#define antPhaseTerms(i,j)              antPhaseBG(beamAnt,i,j)         
#define T_wave                          T_waveBG(beamAnt)               
#define exc_signal                      exc_signalBG(beamAnt)           
#define ant_x                           ant_xBG(beamAnt)               
#define ant_y                           ant_yBG(beamAnt)               
#define ant_z                           ant_zBG(beamAnt)               
#define rampUpMethod                    rampUpMBG(beamAnt)              
#define t_omega                         omega_tBG(beamAnt)              
#define antAngle_zx                     antAngle_zxBG(beamAnt)          
#define antAngle_zy                     antAngle_zyBG(beamAnt)          
#define ant_w0x                         ant_w0xBG(beamAnt)              
#define ant_w0y                         ant_w0yBG(beamAnt)              
#define z2waist                         z2waistBG(beamAnt)              

/*Macros for ABC*/
#define ecoBG(boundaryG)                        boundaryG->eco

#define eco                                     ecoBG(boundaryG)

/*Macros for Mur boundary*/
#define E_Xdir_OLD_G(boundaryG,i,j,k)           boundaryG->E_Xdir_OLD[((i) * (Ny) + j) * (Nz) + k]
#define E_Ydir_OLD_G(boundaryG,i,j,k)           boundaryG->E_Ydir_OLD[((i) * (d_absorb) + j) * (Nz) + k]
#define E_Zdir_OLD_G(boundaryG,i,j,k)           boundaryG->E_Zdir_OLD[((i) * (Ny) + j) * (d_absorb) + k]
#define E_Xdir_OLD_ref_G(boundaryG,i,j,k)       boundaryG->E_Xdir_OLD_ref[((i) * (Ny) + j) * (Nz_ref) + k]
#define E_Ydir_OLD_ref_G(boundaryG,i,j,k)       boundaryG->E_Ydir_OLD_ref[((i) * (d_absorb) + j) * (Nz_ref) + k]
#define E_Zdir_OLD_ref_G(boundaryG,i,j,k)       boundaryG->E_Zdir_OLD_ref[((i) * (Ny) + j) * (d_absorb) + k]

#define E_Xdir_OLD(i,j,k)                       E_Xdir_OLD_G(boundaryG,i,j,k)
#define E_Ydir_OLD(i,j,k)                       E_Ydir_OLD_G(boundaryG,i,j,k)
#define E_Zdir_OLD(i,j,k)                       E_Zdir_OLD_G(boundaryG,i,j,k)
#define E_Xdir_OLD_ref(i,j,k)                   E_Xdir_OLD_ref_G(boundaryG,i,j,k)
#define E_Ydir_OLD_ref(i,j,k)                   E_Ydir_OLD_ref_G(boundaryG,i,j,k)
#define E_Zdir_OLD_ref(i,j,k)                   E_Zdir_OLD_ref_G(boundaryG,i,j,k)

/*Macros for Antenna detector*/
#define DET_ANT_ACCES(antDetect, id, i, j)  \
    ((id == FIELD_01) ? ( antDetect->detAnt_01_fields[ ((i) * 5) + j ] ) : \
     (id == FIELD_02) ? ( antDetect->detAnt_01_fields[ ((i) * 5) + j ] ) : \
     (id == FIELD_03) ? ( antDetect->detAnt_03_fields[ ((i) * 5) + j ] ) : \
     (id == FIELD_04) ? ( antDetect->detAnt_03_fields[ ((i) * 5) + j ] ) : 0 )

#define antDetect01EBG(antDetect,i,j)           antDetect->detAnt_01_fields[ ((i) * 5) + j ]
#define antDetect02EBG(antDetect,i,j)           antDetect->detAnt_02_fields[ ((i) * 5) + j ]
#define antDetect03EBG(antDetect,i,j)           antDetect->detAnt_03_fields[ ((i) * 5) + j ]
#define antDetect04EBG(antDetect,i,j)           antDetect->detAnt_04_fields[ ((i) * 5) + j ]
#define antDetect_1DG(antDetect)                antDetect->antDetect_1D
#define detAnt01zG(antDetect)                   antDetect->detAnt_01_zpos
#define detAnt02zG(antDetect)                   antDetect->detAnt_02_zpos
#define detAnt03zG(antDetect)                   antDetect->detAnt_03_zpos
#define detAnt04zG(antDetect)                   antDetect->detAnt_04_zpos
#define detAnt01yG(antDetect)                   antDetect->detAnt_01_ypos

#define detAnt_01_Fields(i,j)                   antDetect01EBG(antDetect,i,j)
#define detAnt_02_Fields(i,j)                   antDetect02EBG(antDetect,i,j)
#define detAnt_03_Fields(i,j)                   antDetect03EBG(antDetect,i,j)           
#define detAnt_04_Fields(i,j)                   antDetect04EBG(antDetect,i,j)
#define antDetect_1D                            antDetect_1DG(antDetect)
#define detAnt_01_z                             detAnt01zG(antDetect)
#define detAnt_02_z                             detAnt02zG(antDetect)
#define detAnt_03_z                             detAnt03zG(antDetect)
#define detAnt_04_z                             detAnt04zG(antDetect)
#define detAnt_01_y                             detAnt01yG(antDetect)


/*Macros for power struct*/
#define powerDectS(powerValStr)                 powerValStr->pwr_dect
#define powerAbsX1S(powerValStr)                powerValStr->power_abs_x1
#define powerAbsX2S(powerValStr)                powerValStr->power_abs_x2
#define powerAbsY1S(powerValStr)                powerValStr->power_abs_y1
#define powerAbsY2S(powerValStr)                powerValStr->power_abs_y2
#define powerAbsZ1S(powerValStr)                powerValStr->power_abs_z1
#define powerAbsZ2S(powerValStr)                powerValStr->power_abs_z2
#define powerAbsRefS(powerValStr)               powerValStr->power_abs_ref
#define powerPoyX1S(powerValStr)                powerValStr->poynt_x1
#define powerPoyX2S(powerValStr)                powerValStr->poynt_x2
#define powerPoyY1S(powerValStr)                powerValStr->poynt_y1
#define powerPoyY2S(powerValStr)                powerValStr->poynt_y2
#define powerPoyZ1S(powerValStr)                powerValStr->poynt_z1
#define powerPoyZRS(powerValStr)                powerValStr->poynt_z1_ref
#define powerPoyZ2S(powerValStr)                powerValStr->poynt_z2

#define pwr_dect                                powerDectS(powerValStr)
#define power_abs_x1                            powerAbsX1S(powerValStr)
#define power_abs_x2                            powerAbsX2S(powerValStr)                
#define power_abs_y1                            powerAbsY1S(powerValStr)               
#define power_abs_y2                            powerAbsY2S(powerValStr)                
#define power_abs_z1                            powerAbsZ1S(powerValStr)                
#define power_abs_z2                            powerAbsZ2S(powerValStr)                
#define power_abs_ref                           powerAbsRefS(powerValStr)               
#define poynt_x1                                powerPoyX1S(powerValStr)                
#define poynt_x2                                powerPoyX2S(powerValStr)                
#define poynt_y1                                powerPoyY1S(powerValStr)
#define poynt_y2                                powerPoyY2S(powerValStr)                
#define poynt_z1                                powerPoyZ1S(powerValStr)   
#define poynt_z1_ref                            powerPoyZRS(powerValStr) 
#define poynt_z2                                powerPoyZ2S(powerValStr)                

#endif
