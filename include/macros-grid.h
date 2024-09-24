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
#define ne_maxG(gridCfg)        gridCfg->ne_max
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
#define ne_max                  ne_maxG(gridCfg)
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
#define omega_t                         omega_tBG(beamAnt)              
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


#endif
