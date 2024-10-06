#include "auxiliar_module.h"

int set2zero_1D( size_t N_x, double arr_1D[N_x] ){
//{{{

    size_t
        ii;

#pragma omp parallel for default(shared) private(ii)
    for (ii=0 ; ii<N_x ; ++ii) {
        arr_1D[ii] = .0;
    }

    return EXIT_SUCCESS;
} //}}}


int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] ){
//{{{

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii < N_x ; ++ii) {
        for (jj=0 ; jj < N_y ; ++jj) {
            for (kk=0 ; kk < N_z ; ++kk) {
                arr_3D[ii][jj][kk]  = .0;
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}

int setZero2save( gridConfiguration *gridCfg, saveData *saveDCfg){

    size_t
        ii, jj, kk;

    #pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii < Nx/2 ; ++ii) {
        for (jj=0 ; jj < Ny/2 ; ++jj) {
            for (kk=0 ; kk < Nz/2 ; ++kk) {
                data2save(ii,jj,kk)  = .0;
            }
        }
    }

    return EXIT_SUCCESS;

}

int set_densityInAbsorber_v2(   gridConfiguration *gridCfg,
                                systemGrid *G,
                                char absorber[] ) {
//{{{

    double
        x, y, z,
        x0,x1, y0,y1, z0,z1,
        ne_absorb,              // electron density in absorber
        smooth,                 // 1: relatively steep, .2: more smooth
        ne_dist,                // proportional to distance to absorber, where n_e starts to decrease
        scale_fact;

    ne_absorb   = .0;
    smooth      = .5;//.2;
    ne_dist     = round( period/1 );

    x0          = (double)d_absorb + ne_dist;
    x1          = (double)Nx - (d_absorb + ne_dist);
    y0          = (double)d_absorb + ne_dist;
    y1          = (double)Ny - (d_absorb + ne_dist);
    z0          = (double)d_absorb + ne_dist;
    z1          = (double)Nz - (d_absorb + ne_dist);

    // scale to density grid which is only half the size of FDTD-wavefields grid
    // since 2 numerical grid points correspond to one "physical" grid point
    // and we want not to vary the background parameters within one physical
    // grid point
    x0  = round(x0*.5);
    x1  = round(x1*.5);
    y0  = round(y0*.5);
    y1  = round(y1*.5);
    z0  = round(z0*.5);
    z1  = round(z1*.5);

    // the string "absorber" is used to set in which absorber n_e will be modified
    // the comparison is done with the strstr() function, which return the address
    // of the substring if found, NULL otherwise

    // set density in x0 absorber
    if ( strstr(absorber,"x1") ) {
        for ( x=0; x < (Nx/2) ; ++x ) {
            scale_fact  = +.5*(    tanh(smooth*(x-x0)) + 1);        // x0 boundary
            //printf( "x1: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y < (Ny/2) ; ++y )   {
                for ( z=0 ; z < (Nz/2) ; ++z) {
                    n_e((int)x,(int)y,(int)z)  *= scale_fact;
                }
            }
        }
    }
    // set density in x1 absorber
    if ( strstr(absorber,"x2") ) {
        for ( x=0; x<(Nx/2) ; ++x ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(x-x1)) + 1);       // x1 boundary
            //printf( "x2: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(Ny/2) ; ++y )   {
                for (z=0 ; z<(Nz/2) ; ++z) {
                    n_e((int)x,(int)y,(int)z)  *= scale_fact;
                }
            }
        }
    }

    // set density in y0 absorber
    if ( strstr(absorber,"y1") ) {
        for ( y=0; y<(Ny/2) ; ++y ) {
            scale_fact  = +.5*(    tanh(smooth*(y-y0)) + 1);        // y0 boundary
            //printf( "y1: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(Nx/2) ; ++x ) {
                for (z=0 ; z<(Nz/2) ; ++z) {
                    n_e((int)x,(int)y,(int)z)  *= scale_fact;
                }
            }
        }
    }
    // set density in y1 absorber
    if ( strstr(absorber,"y2") ) {
        for ( y=0; y<(Ny/2) ; ++y ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(y-y1)) + 1);       // y1 boundary
            //printf( "y2: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(Nx/2) ; ++x ) {
                for (z=0 ; z<(Nz/2) ; ++z) {
                    n_e((int)x,(int)y,(int)z)  *= scale_fact;
                }
            }
        }
    }

    // set density in z0 absorber
    if ( strstr(absorber,"z1") ) {
        for ( z=0 ; z<(Nz/2) ; ++z) {
            scale_fact  = +.5*(    tanh(smooth*(z-z0)) + 1);        // z0 boundary
            //printf( "z1: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(Nx/2) ; ++x ) {
                for ( y=0; y<(Ny/2) ; ++y ) {
                    n_e((int)x,(int)y,(int)z)  *= scale_fact;
                }
            }
        }
    }
    // set density in z1 absorber
    if ( strstr(absorber,"z2") ) {
        for ( z=0 ; z<(Nz/2) ; ++z) {
            scale_fact  = +.5*(-1.*tanh(smooth*(z-z1)) + 1);       // z1 boundary
            //printf( "z2: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(Nx/2) ; ++x ) {
                for ( y=0; y<(Ny/2) ; ++y ) {
                    n_e((int)x,(int)y,(int)z)  *= scale_fact;
                }
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}