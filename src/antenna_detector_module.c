#include "antenna_detector_module.h"

/*Initialize antenna detector*/
void init_antennaDetect(    antennaDetector *antDetect, 
                            gridConfiguration *gridCfg, 
                            beamAntennaConfiguration *beamAnt){

    initialize_antDetect( antDetect, gridCfg, beamAnt);
    print_antDetect_info( antDetect);


}

void initialize_antDetect(  antennaDetector *antDetect, 
                            gridConfiguration *gridCfg, 
                            beamAntennaConfiguration *beamAnt){

    if( antDetect_1D == 1 ){
        
        ALLOC_2D(antDetect->detAnt_01_fields, Nx, 5, double );
        ALLOC_2D(antDetect->detAnt_02_fields, Nx, 5, double );
        ALLOC_2D(antDetect->detAnt_03_fields, Nx, 5, double );
        ALLOC_2D(antDetect->detAnt_04_fields, Nx, 5, double );

        detAnt_01_y  = ant_y;
        detAnt_01_z  = ant_z + 2;
        detAnt_02_z  = round( ant_z + 2 + 1*5*period ); // steps of 5 cm for 28 GHz = 4.67*period
        detAnt_03_z  = round( ant_z + 2 + 2*5*period );
        detAnt_04_z  = round( ant_z + 2 + 3*5*period );
        // positions have to be even numbers, to ensure fields are accessed correctly
        if ((detAnt_01_y % 2) != 0)  ++detAnt_01_y;
        if ((detAnt_01_z % 2) != 0)  ++detAnt_01_z;
        if ((detAnt_02_z % 2) != 0)  ++detAnt_02_z;
        if ((detAnt_03_z % 2) != 0)  ++detAnt_03_z;
        if ((detAnt_04_z % 2) != 0)  ++detAnt_04_z;
        // issue a warning when detector antenna position is beyond Nz
        if ( detAnt_04_z > ( Nz - d_absorb )) {
            printf( "ERROR: check the detector antenna positions into z direction\n" );
            printf( "       Nz-d_absorb = %d, detAnt_04_zpos = %d\n", 
                    Nz - d_absorb, detAnt_04_z );

        }
    }

    if( antDetect_1D != 1 ){
        printf( "Antenna Detector has not been initialized\n");
    }

}

void print_antDetect_info(antennaDetector *antDetect){

    if( antDetect_1D == 1 ){
    printf( "detector antenna positions: z1 = %d, y1 = %d\n", detAnt_01_z, detAnt_01_y );
    printf( "detector antenna positions: z2 = %d, y1 = %d\n", detAnt_02_z, detAnt_01_y );
    printf( "detector antenna positions: z3 = %d, y1 = %d\n", detAnt_03_z, detAnt_01_y );
    printf( "detector antenna positions: z4 = %d, y1 = %d\n", detAnt_04_z, detAnt_01_y );
    }

}

/*Antenna detect store values*/
void control_antenna_detector( antennaDetector *antDetect, gridConfiguration *gridCfg, systemGrid *G, int t_int ){

    if( antDetect_1D == 1){
        if ( t_int >= ( t_end - 10*period) ) {
            if ( detAnt_01_z < ( Nz - d_absorb )) {
                detAnt_01_storeValues( gridCfg, G, antDetect, detAnt_01_y, detAnt_01_z, t_int );
            }
            if ( detAnt_02_z < ( Nz - d_absorb )) {
                detAnt_02_storeValues( gridCfg, G, antDetect, detAnt_01_y, detAnt_02_z, t_int );
            }
            if ( detAnt_03_z < ( Nz - d_absorb )) {
                detAnt_03_storeValues( gridCfg, G, antDetect, detAnt_01_y, detAnt_03_z, t_int );
            }
            if ( detAnt_04_z < ( Nz - d_absorb )) {
                detAnt_04_storeValues( gridCfg, G, antDetect, detAnt_01_y, detAnt_04_z, t_int );
            }
        }
    }//end if 

}

int detAnt_01_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int ) { 
    //{{{
    size_t
        ii;

    double
        foo;

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

#pragma omp parallel default(shared) private(ii,foo)
#pragma omp for
    for ( ii=2 ; ii <= Nx-2 ; ii+=2 ) {
        // calculate abs(E)
        foo = sqrt(  pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2) );

        // sum of E over time
        // Ex*Ex
        detAnt_01_Fields(ii/2,0) += pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2);
        // Ey*Ey
        detAnt_01_Fields(ii/2,1) += pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2);
        // Ez*Ez
        detAnt_01_Fields(ii/2,2) += pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2);
        // E*E
        detAnt_01_Fields(ii/2,3) += foo*foo;

        // corresponding to an rms(E)-like quantity
        detAnt_01_Fields(ii/2,4) += ( foo * sqrt(1./( (double)(t_int)/(double)(period) + 1e-6 )) );

        //printf( "tt = %d, ii = %d, sum_t(E*E) = %13.5e\n",
        //        tt, ii, detAnt_fields[ii/2][3] );
    }

    return EXIT_SUCCESS;

}//}}}

int detAnt_02_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int ) { 
    //{{{
    size_t
        ii;

    double
        foo;

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

#pragma omp parallel default(shared) private(ii,foo)
#pragma omp for
    for ( ii=2 ; ii <= Nx-2 ; ii+=2 ) {
        // calculate abs(E)
        foo = sqrt(  pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2) );

        // sum of E over time
        // Ex*Ex
        detAnt_02_Fields(ii/2,0) += pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2);
        // Ey*Ey
        detAnt_02_Fields(ii/2,1) += pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2);
        // Ez*Ez
        detAnt_02_Fields(ii/2,2) += pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2);
        // E*E
        detAnt_02_Fields(ii/2,3) += foo*foo;

        // corresponding to an rms(E)-like quantity
        detAnt_02_Fields(ii/2,4) += ( foo * sqrt(1./( (double)(t_int)/(double)(period) + 1e-6 )) );

        //printf( "tt = %d, ii = %d, sum_t(E*E) = %13.5e\n",
        //        tt, ii, detAnt_fields[ii/2][3] );
    }

    return EXIT_SUCCESS;

}//}}}

int detAnt_03_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int ) { 
    //{{{
    size_t
        ii;

    double
        foo;

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

#pragma omp parallel default(shared) private(ii,foo)
#pragma omp for
    for ( ii=2 ; ii <= Nx-2 ; ii+=2 ) {
        // calculate abs(E)
        foo = sqrt(  pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2) );

        // sum of E over time
        // Ex*Ex
        detAnt_03_Fields(ii/2,0) += pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2);
        // Ey*Ey
        detAnt_03_Fields(ii/2,1) += pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2);
        // Ez*Ez
        detAnt_03_Fields(ii/2,2) += pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2);
        // E*E
        detAnt_03_Fields(ii/2,3) += foo*foo;

        // corresponding to an rms(E)-like quantity
        detAnt_03_Fields(ii/2,4) += ( foo * sqrt(1./( (double)(t_int)/(double)(period) + 1e-6 )) );

        //printf( "tt = %d, ii = %d, sum_t(E*E) = %13.5e\n",
        //        tt, ii, detAnt_fields[ii/2][3] );
    }

    return EXIT_SUCCESS;

}//}}}

int detAnt_04_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int ) { 
    //{{{
    size_t
        ii;

    double
        foo;

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

#pragma omp parallel default(shared) private(ii,foo)
#pragma omp for
    for ( ii=2 ; ii <= Nx-2 ; ii+=2 ) {
        // calculate abs(E)
        foo = sqrt(  pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2) );

        // sum of E over time
        // Ex*Ex
        detAnt_04_Fields(ii/2,0) += pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2);
        // Ey*Ey
        detAnt_04_Fields(ii/2,1) += pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2);
        // Ez*Ez
        detAnt_04_Fields(ii/2,2) += pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2);
        // E*E
        detAnt_04_Fields(ii/2,3) += foo*foo;

        // corresponding to an rms(E)-like quantity
        detAnt_04_Fields(ii/2,4) += ( foo * sqrt(1./( (double)(t_int)/(double)(period) + 1e-6 )) );

        //printf( "tt = %d, ii = %d, sum_t(E*E) = %13.5e\n",
        //        tt, ii, detAnt_fields[ii/2][3] );
    }

    return EXIT_SUCCESS;

}//}}}

/*int detAnt1D_storeValues(  gridConfiguration *gridCfg, 
                            systemGrid *G,
                            antennaDetector *antDetect,
                            enum FieldID id,
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int t_int ) { 
    //{{{
    size_t
        ii;

    double
        foo;

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

#pragma omp parallel default(shared) private(ii,foo)
#pragma omp for
    for ( ii=2 ; ii <= Nx-2 ; ii+=2 ) {
        // calculate abs(E)
        foo = sqrt(  pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2)
                    +pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2) );

        // sum of E over time
        // Ex*Ex
        DET_ANT_ACCES(antDetect, id, ii/2, 0) += pow( EB_WAVE(ii+1,detAnt_ypos  ,detAnt_zpos  ), 2);
        // Ey*Ey
        DET_ANT_ACCES(antDetect, id, ii/2, 1) += pow( EB_WAVE(ii  ,detAnt_ypos+1,detAnt_zpos  ), 2);
        // Ez*Ez
        DET_ANT_ACCES(antDetect, id, ii/2, 2) += pow( EB_WAVE(ii  ,detAnt_ypos  ,detAnt_zpos+1), 2);
        // E*E
        DET_ANT_ACCES(antDetect, id, ii/2, 3) += foo*foo;

        // corresponding to an rms(E)-like quantity
        DET_ANT_ACCES(antDetect, id, ii/2, 4) += ( foo * sqrt(1./( (double)(t_int)/(double)(period) + 1e-6 )) );

        //printf( "tt = %d, ii = %d, sum_t(E*E) = %13.5e\n",
        //        tt, ii, detAnt_fields[ii/2][3] );
    }

    return EXIT_SUCCESS;

}//}}}*/