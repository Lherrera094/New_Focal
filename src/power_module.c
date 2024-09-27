#include "power_module.h"

/*Initialize power values*/
void initialize_power( gridConfiguration *gridCfg, powerCalcValues *powerValStr ){
        
    pwr_dect = d_absorb;
    
    printf( "Starting to set all variables to 0...\n" );
    power_abs_x1    = .0;
    power_abs_x2    = .0;
    power_abs_y1    = .0;
    power_abs_y2    = .0;
    power_abs_z1    = .0;
    power_abs_z2    = .0;
    power_abs_ref   = 1e-7;
    poynt_x1       = .0;
    poynt_x2       = .0;
    poynt_y1       = .0;
    poynt_y2       = .0;
    poynt_z1       = .0;
    poynt_z1_ref   = .0;
    poynt_z2       = .0;
    printf( "...done setting all variables to 0\n" );  
}

/*Functions to compute Poynting vector*/
void control_power(     gridConfiguration *gridCfg, 
                        systemGrid *G, 
                        powerCalcValues *powerValStr,
                        saveData *saveDCfg, 
                        beamAntennaConfiguration *beamAnt,
                        int t_int ){

        compute_power( gridCfg, G, powerValStr, t_int );
        compute_timetraces( gridCfg, saveDCfg, beamAnt, powerValStr, t_int );
}


void compute_power( gridConfiguration *gridCfg, systemGrid *G, powerCalcValues *powerValStr, int t_int ){
    //IQ detector for power detection
    if ( t_int >= 20*period ) {
        // z1-plane and z2-plane
        poynt_z1_ref    = calc_poynt_4( gridCfg, G, powerValStr, "ref_z1" );
        poynt_z1        = calc_poynt_4( gridCfg, G, powerValStr, "z1"     );
        poynt_z2        = calc_poynt_4( gridCfg, G, powerValStr, "z2"     );
        // x1-plane and x2-plane
        poynt_x1        = calc_poynt_4( gridCfg, G, powerValStr, "x1"     );
        poynt_x2        = calc_poynt_4( gridCfg, G, powerValStr, "x2"     );
        // y1-plane and y2-plane
        poynt_y1        = calc_poynt_4( gridCfg, G, powerValStr, "y1"     );
        poynt_y2        = calc_poynt_4( gridCfg, G, powerValStr, "y2"     );

            
//     printf( "t = %d, power_abs_ref = %13.5e, power_abs_z1 = %13.5e, power_abs_z2 = %13.5e, poynt_z1 = %13.5e, poynt_z2 = %13.5e\n",
//              t_int, power_abs_ref, power_abs_z1, power_abs_z2, poynt_z1, poynt_z2 );

        power_abs_ref   = .99*power_abs_ref + .01*poynt_z1_ref;
        power_abs_z1    = .99*power_abs_z1  + .01*poynt_z1;
        power_abs_z2    = .99*power_abs_z2  + .01*poynt_z2;
        power_abs_x1    = .99*power_abs_x1  + .01*poynt_x1;
        power_abs_x2    = .99*power_abs_x2  + .01*poynt_x2;
        power_abs_y1    = .99*power_abs_y1  + .01*poynt_y1;
        power_abs_y2    = .99*power_abs_y2  + .01*poynt_y2;
    }

}

double calc_poynt_4(    gridConfiguration *gridCfg,
                        systemGrid *G,
                        powerCalcValues *powerValStr,
                        char absorber[] ) {
//{{{

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii = pwr_dect ; ii <= (Nx-pwr_dect-2) ; ii+=2) {
            for (jj = pwr_dect ; jj <= (Ny-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref(ii+1,jj  ,pwr_dect  )
                          *EB_WAVE_ref(ii+1,jj  ,pwr_dect+1)
                          -EB_WAVE_ref(ii  ,jj+1,pwr_dect  )
                          *EB_WAVE_ref(ii  ,jj+1,pwr_dect+1) );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii = pwr_dect ; ii <= (Nx-pwr_dect-2) ; ii+=2) {
            for (jj = pwr_dect ; jj <= (Ny-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE(ii+1,jj  ,pwr_dect  ) - EB_WAVE_ref(ii+1,jj  ,pwr_dect  ) )
                          *( EB_WAVE(ii+1,jj  ,pwr_dect+1) - EB_WAVE_ref(ii+1,jj  ,pwr_dect+1) )
                          -( EB_WAVE(ii  ,jj+1,pwr_dect  ) - EB_WAVE_ref(ii  ,jj+1,pwr_dect  ) ) 
                          *( EB_WAVE(ii  ,jj+1,pwr_dect+1) - EB_WAVE_ref(ii  ,jj+1,pwr_dect+1) ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii = pwr_dect ; ii <= (Nx-pwr_dect-2) ; ii+=2) {
            for (jj = pwr_dect ; jj <= (Ny-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE(ii+1,jj  ,Nz-pwr_dect  )
                          *EB_WAVE(ii+1,jj  ,Nz-pwr_dect-1)
                          -EB_WAVE(ii  ,jj+1,Nz-pwr_dect  )
                          *EB_WAVE(ii  ,jj+1,Nz-pwr_dect-1) );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj = pwr_dect ; jj <= (Ny-pwr_dect-2) ; jj+=2) {
            for (kk = pwr_dect ; kk <= (Nz-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE(pwr_dect  ,jj+1,kk  )
                          *EB_WAVE(pwr_dect+1,jj+1,kk  )
                          -EB_WAVE(pwr_dect  ,jj  ,kk+1)
                          *EB_WAVE(pwr_dect+1,jj  ,kk+1) );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj = pwr_dect ; jj <= (Ny-pwr_dect-2) ; jj+=2) {
            for (kk = pwr_dect ; kk <= (Nz-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE(Nx-pwr_dect  ,jj+1,kk  )
                          *EB_WAVE(Nx-pwr_dect-1,jj+1,kk  )
                          -EB_WAVE(Nx-pwr_dect  ,jj  ,kk+1)
                          *EB_WAVE(Nx-pwr_dect-1,jj  ,kk+1) );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii = pwr_dect ; ii <= (Nx-pwr_dect-2) ; ii+=2) {
            for (kk = pwr_dect ; kk <= (Nz-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE(ii  ,pwr_dect  ,kk+1)
                          *EB_WAVE(ii  ,pwr_dect+1,kk+1)
                          -EB_WAVE(ii+1,pwr_dect  ,kk  )
                          *EB_WAVE(ii+1,pwr_dect+1,kk  ) );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii = pwr_dect ; ii <= (Nx-pwr_dect-2) ; ii+=2) {
            for (kk = pwr_dect ; kk <= (Nz-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE(ii  ,Ny-pwr_dect  ,kk+1)
                          *EB_WAVE(ii  ,Ny-pwr_dect-1,kk+1)
                          -EB_WAVE(ii+1,Ny-pwr_dect  ,kk  )
                          *EB_WAVE(ii+1,Ny-pwr_dect-1,kk  ) );
            }
        }
    }
    
    return fabs(poynt);
} //}}}

/*Timetraces computing functions*/
void compute_timetraces(        gridConfiguration *gridCfg, 
                                saveData *saveDCfg, 
                                beamAntennaConfiguration *beamAnt,
                                powerCalcValues *powerValStr,
                                int t_int ){

    if ( (t_int % (int)period) == 4 )  {
            printf( "status: number of oscillation periods: %d (t_int= %d) \n",T_wave,t_int);
            printf( "        Poynting-power: z1 = %13.6e, z2 = %13.6e, x1 = %13.6e, x2 = %13.6e, y1 = %13.6e, y2 = %13.6e, (z1+z2+x1+x2+y1+y2)/z1_ref = %13.6e %%\n",
                    power_abs_z1/power_abs_ref, 
                    power_abs_z2/power_abs_ref,
                    power_abs_x1/power_abs_ref, 
                    power_abs_x2/power_abs_ref,
                    power_abs_y1/power_abs_ref, 
                    power_abs_y2/power_abs_ref,
                    (power_abs_x1+power_abs_x2 + power_abs_y1+power_abs_y2 + power_abs_z1+power_abs_z2)/power_abs_ref * 100.
                    );
            timetraces(T_wave,0)   = (double)t_int;
            timetraces(T_wave,1)   = (double)T_wave;
            timetraces(T_wave,2)   = power_abs_z1/power_abs_ref;
            timetraces(T_wave,3)   = power_abs_z2/power_abs_ref;
            timetraces(T_wave,4)   = power_abs_x1/power_abs_ref;
            timetraces(T_wave,5)   = power_abs_x2/power_abs_ref;
            timetraces(T_wave,6)   = power_abs_y1/power_abs_ref;
            timetraces(T_wave,7)   = power_abs_y2/power_abs_ref;

        }

}