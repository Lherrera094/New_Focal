#include "boundary_module.h"

/*Initialize boundary conditions*/
void init_boundary(gridConfiguration *gridCfg, boundaryGrid *boundaryG){

    if(boundary_sel == 1){

        eco = 10./(double)(period);

    }
    else if (boundary_sel == 2){

        ALLOC_3D(boundaryG->E_Xdir_OLD, d_absorb, Ny, Nz, double);
        ALLOC_3D(boundaryG->E_Ydir_OLD, Nx, d_absorb, Nz, double);
        ALLOC_3D(boundaryG->E_Zdir_OLD, Nx, Ny, d_absorb, double);

        ALLOC_3D(boundaryG->E_Xdir_OLD_ref, d_absorb, Ny, Nz_ref, double);
        ALLOC_3D(boundaryG->E_Ydir_OLD_ref, Nx, d_absorb, Nz_ref, double);
        ALLOC_3D(boundaryG->E_Zdir_OLD_ref, Nx, Ny, d_absorb,     double);
    }
    else if(boundary_sel == 3){

        ALLOC_3D(boundaryG->DH_WAVE, Nx, Ny, Nz, double);
        ALLOC_3D(boundaryG->DH_WAVE_ref, Nx, Ny, Nz_ref, double);

        ALLOC_1D(boundaryG->F1x, Nx/2, double);
        ALLOC_1D(boundaryG->F1y, Ny/2, double);
        ALLOC_1D(boundaryG->F1z, Nz/2, double);
        ALLOC_1D(boundaryG->F2x, Nx/2, double);
        ALLOC_1D(boundaryG->F2y, Ny/2, double);
        ALLOC_1D(boundaryG->F2z, Nz/2, double);
        ALLOC_1D(boundaryG->Cx, Nx/2, double);
        ALLOC_1D(boundaryG->Cy, Ny/2, double);
        ALLOC_1D(boundaryG->Cz, Nz/2, double);

        /*UPML parameters for reference*/
        ALLOC_1D(boundaryG->F1zr, Nz_ref/2, double);
        ALLOC_1D(boundaryG->F2zr, Nz_ref/2, double);
        ALLOC_1D(boundaryG->Czr, Nz_ref/2, double);
        
        init_UPML_parameters(   gridCfg, boundaryG);

    }
   
}

/*Apply boundary on time evolution*/
void advance_boundary(gridConfiguration *gridCfg, systemGrid *G, boundaryGrid *boundaryG){

    if(boundary_sel == 1){

        apply_absorber( gridCfg, G, boundaryG);
        apply_absorber_ref(gridCfg, G, boundaryG);

    }
    else if (boundary_sel == 2){
        
        abc_Mur_1st(gridCfg,"x1x2y1y2z1z2", G, boundaryG);
        abc_Mur_saveOldE_xdir(gridCfg, G, boundaryG);
        abc_Mur_saveOldE_ydir(gridCfg, G, boundaryG);
        abc_Mur_saveOldE_zdir(gridCfg, G, boundaryG);
        
        abc_Mur_1st_ref(gridCfg, G, boundaryG);
        abc_Mur_saveOldE_ref_xdir(gridCfg, G, boundaryG);
        abc_Mur_saveOldE_ref_ydir(gridCfg, G, boundaryG);
        abc_Mur_saveOldE_ref_zdir(gridCfg, G, boundaryG);
        
    }

}

/*Section for simple Absorbing boundary conditions*/
int apply_absorber( gridConfiguration *gridCfg, 
                    systemGrid *G, 
                    boundaryGrid *boundaryG ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb-2 ; kk+=2) {

                damp = ((double)kk-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE(ii  ,jj  ,kk+1) *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)

    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk = (Nz - d_absorb) ; kk < Nz-2 ; kk+=2) {      //Nz-d_absorb-2 ???

                damp = ((double)kk-((double)Nz-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < d_absorb-2 ; ii+=2){
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                damp = ((double)ii-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii = (Nx - d_absorb) ; ii < Nx-2 ; ii+=2){ //Nx-d_absorb-2 ???
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {    
                damp = ((double)ii-((double)Nx-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                damp = ((double)jj-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj = (Ny - d_absorb) ; jj < Ny-2 ; jj+=2) { //Ny-d_absorb-2 ???
            for (kk=2 ; kk < Nz-2 ; kk+=2) {  
                damp = ((double)jj-((double)Ny-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}

int apply_absorber_ref( gridConfiguration *gridCfg, 
                        systemGrid *G, 
                        boundaryGrid *boundaryG ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb-2 ; kk+=2) {
                damp = ((double)kk-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE_ref(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj  ,kk+1) *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk = (Nz_ref - d_absorb) ; kk < Nz_ref-2 ; kk+=2) {      //Nz-d_absorb-2 ???
                damp = ((double)kk-((double)Nz_ref-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE_ref(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj < Ny-2 ; jj+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            for (ii=2 ; ii < d_absorb-2 ; ii+=2) {
                damp = ((double)ii-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE_ref(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj < Ny-2 ; jj+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {  
            for (ii = (Nx - d_absorb) ; ii < Nx-2 ; ii+=2) {    //Nx-d_absorb-2 ???
                damp = ((double)ii-((double)Nx-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE_ref(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            for (jj=2 ; jj < d_absorb-2 ; jj+=2) {
                damp = ((double)jj-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE_ref(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            for (jj = (Ny - d_absorb) ; jj < Ny-2 ; jj+=2) {  //Ny-d_absorb-2 ???
                damp = ((double)jj-((double)Ny-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE_ref(ii+1,jj  ,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj+1,kk  ) *= damp;
                EB_WAVE_ref(ii  ,jj  ,kk+1) *= damp;
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}

/*Section for Mur boundary implementation*/
int abc_Mur_saveOldE_xdir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj < Ny-2 ; jj+=2) {
        for (kk=2 ; kk < Nz-2 ; kk+=2) {
            // store values at x=0 and x=1
            // Ex: odd-even-even
            E_Xdir_OLD(0+1,jj  ,kk  ) = EB_WAVE(0+offset+1,jj  ,kk  );
            E_Xdir_OLD(2+1,jj  ,kk  ) = EB_WAVE(2+offset+1,jj  ,kk  );
            // Ey: even-odd-even
            E_Xdir_OLD(0  ,jj+1,kk  ) = EB_WAVE(0+offset  ,jj+1,kk  );
            E_Xdir_OLD(2  ,jj+1,kk  ) = EB_WAVE(2+offset  ,jj+1,kk  );
            // Ez: even-even-odd
            E_Xdir_OLD(0  ,jj  ,kk+1) = EB_WAVE(0+offset  ,jj  ,kk+1);
            E_Xdir_OLD(2  ,jj  ,kk+1) = EB_WAVE(2+offset  ,jj  ,kk+1);

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_Xdir_OLD(4+1,jj  ,kk  ) = EB_WAVE(Nx-4-offset+1,jj  ,kk  );
            E_Xdir_OLD(6+1,jj  ,kk  ) = EB_WAVE(Nx-2-offset+1,jj  ,kk  );
            // Ey: even-odd-even
            E_Xdir_OLD(4  ,jj+1,kk  ) = EB_WAVE(Nx-4-offset  ,jj+1,kk  );
            E_Xdir_OLD(6  ,jj+1,kk  ) = EB_WAVE(Nx-2-offset  ,jj+1,kk  );
            // Ez: even-even-odd
            E_Xdir_OLD(4  ,jj  ,kk+1) = EB_WAVE(Nx-4-offset  ,jj  ,kk+1);
            E_Xdir_OLD(6  ,jj  ,kk+1) = EB_WAVE(Nx-2-offset  ,jj  ,kk+1);
        }
    }
 
    return EXIT_SUCCESS;
}//}}}

int abc_Mur_saveOldE_ydir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (kk=2 ; kk < Nz-2 ; kk+=2) {
            // store values at y=0 and y=1
            // Ex: odd-even-even
            E_Ydir_OLD(ii+1,0  ,kk  ) = EB_WAVE(ii+1,0+offset  ,kk  );
            E_Ydir_OLD(ii+1,2  ,kk  ) = EB_WAVE(ii+1,2+offset  ,kk  );
            // Ey: even-odd-even
            E_Ydir_OLD(ii  ,0+1,kk  ) = EB_WAVE(ii  ,0+offset+1,kk  );
            E_Ydir_OLD(ii  ,2+1,kk  ) = EB_WAVE(ii  ,2+offset+1,kk  );
            // Ez: even-even-odd
            E_Ydir_OLD(ii  ,0  ,kk+1) = EB_WAVE(ii  ,0+offset  ,kk+1);
            E_Ydir_OLD(ii  ,2  ,kk+1) = EB_WAVE(ii  ,2+offset  ,kk+1);

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_Ydir_OLD(ii+1,4  ,kk  ) = EB_WAVE(ii+1,Ny-4-offset  ,kk  );
            E_Ydir_OLD(ii+1,6  ,kk  ) = EB_WAVE(ii+1,Ny-2-offset  ,kk  );
            // Ey: even-odd-even
            E_Ydir_OLD(ii  ,4+1,kk  ) = EB_WAVE(ii  ,Ny-4-offset+1,kk  );
            E_Ydir_OLD(ii  ,6+1,kk  ) = EB_WAVE(ii  ,Ny-2-offset+1,kk  );
            // Ez: even-even-odd
            E_Ydir_OLD(ii  ,4  ,kk+1) = EB_WAVE(ii  ,Ny-4-offset  ,kk+1);
            E_Ydir_OLD(ii  ,6  ,kk+1) = EB_WAVE(ii  ,Ny-2-offset  ,kk+1);
        }
    }
 
    return EXIT_SUCCESS;

}//}}}

int abc_Mur_saveOldE_zdir(  gridConfiguration *gridCfg, 
                            systemGrid *G, 
                            boundaryGrid *boundaryG ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            // store values at z=0 and z=1
            // Ex: odd-even-even
            E_Zdir_OLD(ii+1,jj  ,0  ) = EB_WAVE(ii+1,jj  ,0+offset  );
            E_Zdir_OLD(ii+1,jj  ,2  ) = EB_WAVE(ii+1,jj  ,2+offset  );
            // Ey: even-odd-even
            E_Zdir_OLD(ii  ,jj+1,0  ) = EB_WAVE(ii  ,jj+1,0+offset  );
            E_Zdir_OLD(ii  ,jj+1,2  ) = EB_WAVE(ii  ,jj+1,2+offset  );
            // Ez: even-even-odd
            E_Zdir_OLD(ii  ,jj  ,0+1) = EB_WAVE(ii  ,jj  ,0+offset+1);
            E_Zdir_OLD(ii  ,jj  ,2+1) = EB_WAVE(ii  ,jj  ,2+offset+1);

            // store values at z=Nz-1 and z=Nz-2
            // Ex: odd-even-even
            E_Zdir_OLD(ii+1,jj  ,4  ) = EB_WAVE(ii+1,jj  ,Nz-4-offset  );
            E_Zdir_OLD(ii+1,jj  ,6  ) = EB_WAVE(ii+1,jj  ,Nz-2-offset  );
            // Ey: even-odd-even
            E_Zdir_OLD(ii  ,jj+1,4  ) = EB_WAVE(ii  ,jj+1,Nz-4-offset  );
            E_Zdir_OLD(ii  ,jj+1,6  ) = EB_WAVE(ii  ,jj+1,Nz-2-offset  );
            // Ez: even-even-odd
            E_Zdir_OLD(ii  ,jj  ,4+1) = EB_WAVE(ii  ,jj  ,Nz-4-offset+1);
            E_Zdir_OLD(ii  ,jj  ,6+1) = EB_WAVE(ii  ,jj  ,Nz-2-offset+1);

        }
    }
 
    return EXIT_SUCCESS;

}//}}}

/*Section for Mur boundary implementation*/
int abc_Mur_saveOldE_ref_xdir(  gridConfiguration *gridCfg, 
                                systemGrid *G, 
                                boundaryGrid *boundaryG ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj < Ny-2 ; jj+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            // store values at x=0 and x=1
            // Ex: odd-even-even
            E_Xdir_OLD_ref(0+1,jj  ,kk  ) = EB_WAVE_ref(0+offset+1,jj  ,kk  );
            E_Xdir_OLD_ref(2+1,jj  ,kk  ) = EB_WAVE_ref(2+offset+1,jj  ,kk  );
            // Ey: even-odd-even
            E_Xdir_OLD_ref(0  ,jj+1,kk  ) = EB_WAVE_ref(0+offset  ,jj+1,kk  );
            E_Xdir_OLD_ref(2  ,jj+1,kk  ) = EB_WAVE_ref(2+offset  ,jj+1,kk  );
            // Ez: even-even-odd
            E_Xdir_OLD_ref(0  ,jj  ,kk+1) = EB_WAVE_ref(0+offset  ,jj  ,kk+1);
            E_Xdir_OLD_ref(2  ,jj  ,kk+1) = EB_WAVE_ref(2+offset  ,jj  ,kk+1);

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_Xdir_OLD_ref(4+1,jj  ,kk  ) = EB_WAVE_ref(Nx-4-offset+1,jj  ,kk  );
            E_Xdir_OLD_ref(6+1,jj  ,kk  ) = EB_WAVE_ref(Nx-2-offset+1,jj  ,kk  );
            // Ey: even-odd-even
            E_Xdir_OLD_ref(4  ,jj+1,kk  ) = EB_WAVE_ref(Nx-4-offset  ,jj+1,kk  );
            E_Xdir_OLD_ref(6  ,jj+1,kk  ) = EB_WAVE_ref(Nx-2-offset  ,jj+1,kk  );
            // Ez: even-even-odd
            E_Xdir_OLD_ref(4  ,jj  ,kk+1) = EB_WAVE_ref(Nx-4-offset  ,jj  ,kk+1);
            E_Xdir_OLD_ref(6  ,jj  ,kk+1) = EB_WAVE_ref(Nx-2-offset  ,jj  ,kk+1);
        }
    }
 
    return EXIT_SUCCESS;

}//}}}

int abc_Mur_saveOldE_ref_ydir(  gridConfiguration *gridCfg, 
                                systemGrid *G, 
                                boundaryGrid *boundaryG ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            // store values at y=0 and y=1
            // Ex: odd-even-even
            E_Ydir_OLD_ref(ii+1,0  ,kk  ) = EB_WAVE_ref(ii+1,0+offset  ,kk  );
            E_Ydir_OLD_ref(ii+1,2  ,kk  ) = EB_WAVE_ref(ii+1,2+offset  ,kk  );
            // Ey: even-odd-even
            E_Ydir_OLD_ref(ii  ,0+1,kk  ) = EB_WAVE_ref(ii  ,0+offset+1,kk  );
            E_Ydir_OLD_ref(ii  ,2+1,kk  ) = EB_WAVE_ref(ii  ,2+offset+1,kk  );
            // Ez: even-even-odd
            E_Ydir_OLD_ref(ii  ,0  ,kk+1) = EB_WAVE_ref(ii  ,0+offset  ,kk+1);
            E_Ydir_OLD_ref(ii  ,2  ,kk+1) = EB_WAVE_ref(ii  ,2+offset  ,kk+1);

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_Ydir_OLD_ref(ii+1,4  ,kk  ) = EB_WAVE_ref(ii+1,Ny-4-offset  ,kk  );
            E_Ydir_OLD_ref(ii+1,6  ,kk  ) = EB_WAVE_ref(ii+1,Ny-2-offset  ,kk  );
            // Ey: even-odd-even
            E_Ydir_OLD_ref(ii  ,4+1,kk  ) = EB_WAVE_ref(ii  ,Ny-4-offset+1,kk  );
            E_Ydir_OLD_ref(ii  ,6+1,kk  ) = EB_WAVE_ref(ii  ,Ny-2-offset+1,kk  );
            // Ez: even-even-odd
            E_Ydir_OLD_ref(ii  ,4  ,kk+1) = EB_WAVE_ref(ii  ,Ny-4-offset  ,kk+1);
            E_Ydir_OLD_ref(ii  ,6  ,kk+1) = EB_WAVE_ref(ii  ,Ny-2-offset  ,kk+1);
        }
    }
 
    return EXIT_SUCCESS;

}//}}}

int abc_Mur_saveOldE_ref_zdir(  gridConfiguration *gridCfg, 
                                systemGrid *G, 
                                boundaryGrid *boundaryG ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            // store values at z=0 and z=1
            // Ex: odd-even-even
            E_Zdir_OLD_ref(ii+1,jj  ,0  ) = EB_WAVE_ref(ii+1,jj  ,0+offset  );
            E_Zdir_OLD_ref(ii+1,jj  ,2  ) = EB_WAVE_ref(ii+1,jj  ,2+offset  );
            // Ey: even-odd-even
            E_Zdir_OLD_ref(ii  ,jj+1,0  ) = EB_WAVE_ref(ii  ,jj+1,0+offset  );
            E_Zdir_OLD_ref(ii  ,jj+1,2  ) = EB_WAVE_ref(ii  ,jj+1,2+offset  );
            // Ez: even-even-odd
            E_Zdir_OLD_ref(ii  ,jj  ,0+1) = EB_WAVE_ref(ii  ,jj  ,0+offset+1);
            E_Zdir_OLD_ref(ii  ,jj  ,2+1) = EB_WAVE_ref(ii  ,jj  ,2+offset+1);

            // store values at z=Nz-1 and z=Nz-2
            // Ex: odd-even-even
            E_Zdir_OLD_ref(ii+1,jj  ,4  ) = EB_WAVE_ref(ii+1,jj  ,Nz_ref-4-offset  );
            E_Zdir_OLD_ref(ii+1,jj  ,6  ) = EB_WAVE_ref(ii+1,jj  ,Nz_ref-2-offset  );
            // Ey: even-odd-even
            E_Zdir_OLD_ref(ii  ,jj+1,4  ) = EB_WAVE_ref(ii  ,jj+1,Nz_ref-4-offset  );
            E_Zdir_OLD_ref(ii  ,jj+1,6  ) = EB_WAVE_ref(ii  ,jj+1,Nz_ref-2-offset  );
            // Ez: even-even-odd
            E_Zdir_OLD_ref(ii  ,jj  ,4+1) = EB_WAVE_ref(ii  ,jj  ,Nz_ref-4-offset+1);
            E_Zdir_OLD_ref(ii  ,jj  ,6+1) = EB_WAVE_ref(ii  ,jj  ,Nz_ref-2-offset+1);

        }
    }
 
    return EXIT_SUCCESS;

}//}}}

int abc_Mur_1st(    gridConfiguration *gridCfg, 
                    char absorber[],
                    systemGrid *G, 
                    boundaryGrid *boundaryG) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (dt - dx)/(dt + dx);
    offset  = 2;

    // the string "absorber" is used to set which absorber is treated
    // the comparison is done with the strstr() function, which return the address
    // of the substring if found, NULL otherwise
    // NOTE: "if (strstr(absorber,"x1))" should be sufficient
    //       "if (strstr(absorber,"x1) != NULL)" should be equivalent

    // absorber into x-direction
    if ( strstr(absorber,"x1") ) {      
        //printf("abs_Mur_1st_v2: x1\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                // absorber at x=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE(offset+0+1,jj  ,kk  ) = E_Xdir_OLD(2+1,jj  ,kk  )
                                              + cnst * ( EB_WAVE(offset+2+1,jj  ,kk  ) 
                                                     -E_Xdir_OLD(0+1       ,jj  ,kk  ) );
                // Ey: even-odd-even
                EB_WAVE(offset+0  ,jj+1,kk  ) = E_Xdir_OLD(2  ,jj+1,kk  )
                                              + cnst * ( EB_WAVE(offset+2  ,jj+1,kk  ) 
                                                     -E_Xdir_OLD(0         ,jj+1,kk  ) );
                // Ez: even-even-odd
                EB_WAVE(offset+0  ,jj  ,kk+1) = E_Xdir_OLD(2  ,jj  ,kk+1)
                                              + cnst * ( EB_WAVE(offset+2  ,jj  ,kk+1) 
                                                     -E_Xdir_OLD(0         ,jj  ,kk+1) );
            }
        }
    }
    if ( strstr(absorber,"x2") ) {
        //printf("abs_Mur_1st_v2: x2\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                // absorber at x=Nx grid boundary
                // Ex: odd-even-even
                EB_WAVE(Nx-2-offset+1,jj  ,kk  ) = E_Xdir_OLD(4+1,jj  ,kk  )
                                                 + cnst * ( EB_WAVE(Nx-4-offset+1,jj  ,kk  ) 
                                                        -E_Xdir_OLD(6+1          ,jj  ,kk  ) );
                // Ey: even-odd-even
                EB_WAVE(Nx-2-offset  ,jj+1,kk  ) = E_Xdir_OLD(4  ,jj+1,kk  )
                                                 + cnst * ( EB_WAVE(Nx-4-offset  ,jj+1,kk  ) 
                                                        -E_Xdir_OLD(6            ,jj+1,kk  ) );
                // Ez: even-even-odd
                EB_WAVE(Nx-2-offset  ,jj  ,kk+1) = E_Xdir_OLD(4  ,jj  ,kk+1)
                                                 + cnst * ( EB_WAVE(Nx-4-offset  ,jj  ,kk+1) 
                                                        -E_Xdir_OLD(6            ,jj  ,kk+1) );
            }
        }
    }

    // absorber into y-direction
    if ( strstr(absorber,"y1") ) {
        //printf("abs_Mur_1st_v2: y1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii < Nx-2 ; ii+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                // absorber at y=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE(ii+1,offset+0  ,kk  ) = E_Ydir_OLD(ii+1,2  ,kk  )
                                              + cnst * ( EB_WAVE(ii+1,offset+2  ,kk  )
                                                     -E_Ydir_OLD(ii+1,0         ,kk  ) );
                // Ey: even-odd-even
                EB_WAVE(ii  ,offset+0+1,kk  ) = E_Ydir_OLD(ii  ,2+1,kk  )
                                              + cnst * ( EB_WAVE(ii  ,offset+2+1,kk  )
                                                     -E_Ydir_OLD(ii  ,0+1       ,kk  ) );
                // Ez: even-even-odd
                EB_WAVE(ii  ,offset+0  ,kk+1) = E_Ydir_OLD(ii  ,2  ,kk+1)
                                              + cnst * ( EB_WAVE(ii  ,offset+2  ,kk+1)
                                                     -E_Ydir_OLD(ii  ,0         ,kk+1) );
            }
        }
    }
    if ( strstr(absorber,"y2") ) {
        //printf("abs_Mur_1st_v2: y2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii < Nx-2 ; ii+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                // absorber at y=Ny grid boundary
                // Ex: odd-even-even
                EB_WAVE(ii+1,Ny-2-offset  ,kk  ) = E_Ydir_OLD(ii+1,4  ,kk  )
                                                 + cnst * ( EB_WAVE(ii+1,Ny-4-offset  ,kk  )
                                                        -E_Ydir_OLD(ii+1,6            ,kk  ) );
                // Ey: even-odd-even
                EB_WAVE(ii  ,Ny-2-offset+1,kk  ) = E_Ydir_OLD(ii  ,4+1,kk  )
                                                 + cnst * ( EB_WAVE(ii  ,Ny-4-offset+1,kk  )
                                                        -E_Ydir_OLD(ii  ,6+1          ,kk  ) );
                // Ez: even-even-odd
                EB_WAVE(ii  ,Ny-2-offset  ,kk+1) = E_Ydir_OLD(ii  ,4  ,kk+1)
                                                 + cnst * ( EB_WAVE(ii  ,Ny-4-offset  ,kk+1)
                                                        -E_Ydir_OLD(ii  ,6            ,kk+1) );
            }
        }
    }

    // absorber into z-direction
    if ( strstr(absorber,"z1") ) {
        //printf("abs_Mur_1st_v2: z1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii < Nx-2 ; ii+=2) {
            for (jj=2 ; jj < Ny-2 ; jj+=2) {
                // absorber at z=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE(ii+1,jj  ,offset+0)   = E_Zdir_OLD(ii+1,jj  ,2  )
                                              + cnst * ( EB_WAVE(ii+1,jj  ,offset+2  )
                                                     -E_Zdir_OLD(ii+1,jj  ,0         ) );
                // Ey: even-odd-even
                EB_WAVE(ii  ,jj+1,offset+0)   = E_Zdir_OLD(ii  ,jj+1,2  )
                                              + cnst * ( EB_WAVE(ii  ,jj+1,offset+2  )
                                                     -E_Zdir_OLD(ii  ,jj+1,0         ) );
                // Ez: even-even-odd
                EB_WAVE(ii  ,jj  ,offset+0+1) = E_Zdir_OLD(ii  ,jj  ,2+1)
                                              + cnst * ( EB_WAVE(ii  ,jj  ,offset+2+1)
                                                     -E_Zdir_OLD(ii  ,jj  ,0+1       ) );
            }
        }
    }
    if ( strstr(absorber,"z2") ) {
        //printf("abs_Mur_1st_v2: z2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii < Nx-2 ; ii+=2) {
            for (jj=2 ; jj < Ny-2 ; jj+=2) {
                // absorber at z=Nz grid boundary
                // Ex: odd-even-even
                EB_WAVE(ii+1,jj  ,Nz-2-offset  )    = E_Zdir_OLD(ii+1,jj  ,4  )
                                                    + cnst * ( EB_WAVE(ii+1,jj  ,Nz-4-offset  )
                                                           -E_Zdir_OLD(ii+1,jj  ,6            ) );
                // Ey: even-odd-even
                EB_WAVE(ii  ,jj+1,Nz-2-offset  )    = E_Zdir_OLD(ii  ,jj+1,4  )
                                                    + cnst * ( EB_WAVE(ii  ,jj+1,Nz-4-offset  )
                                                           -E_Zdir_OLD(ii  ,jj+1,6            ) );
                // Ez: even-even-odd
                EB_WAVE(ii  ,jj  ,Nz-2-offset+1)    = E_Zdir_OLD(ii  ,jj  ,4+1)
                                                    + cnst * ( EB_WAVE(ii  ,jj  ,Nz-4-offset+1)
                                                           -E_Zdir_OLD(ii  ,jj  ,6+1          ) );
            }
        }
    }

    return EXIT_SUCCESS;

} //}}}

int abc_Mur_1st_ref(    gridConfiguration *gridCfg, 
                        systemGrid *G, 
                        boundaryGrid *boundaryG) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (dt - dx)/(dt + dx);
    offset  = 2;

    // the string "absorber" is used to set which absorber is treated
    // the comparison is done with the strstr() function, which return the address
    // of the substring if found, NULL otherwise
    // NOTE: "if (strstr(absorber,"x1))" should be sufficient
    //       "if (strstr(absorber,"x1) != NULL)" should be equivalent

    // absorber into x-direction
         
    //printf("abs_Mur_1st_v2: x1\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj < Ny-2 ; jj+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            // absorber at x=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE_ref(offset+0+1,jj  ,kk  ) = E_Xdir_OLD_ref(2+1,jj  ,kk  )
                                              + cnst * ( EB_WAVE_ref(offset+2+1,jj  ,kk  ) 
                                                     -E_Xdir_OLD_ref(0+1       ,jj  ,kk  ) );
            // Ey: even-odd-even
            EB_WAVE_ref(offset+0  ,jj+1,kk  ) = E_Xdir_OLD_ref(2  ,jj+1,kk  )
                                              + cnst * ( EB_WAVE_ref(offset+2  ,jj+1,kk  ) 
                                                     -E_Xdir_OLD_ref(0         ,jj+1,kk  ) );
            // Ez: even-even-odd
            EB_WAVE_ref(offset+0  ,jj  ,kk+1) = E_Xdir_OLD_ref(2  ,jj  ,kk+1)
                                              + cnst * ( EB_WAVE_ref(offset+2  ,jj  ,kk+1) 
                                                     -E_Xdir_OLD_ref(0         ,jj  ,kk+1) );
        }
    }
    
    //printf("abs_Mur_1st_v2: x2\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj < Ny-2 ; jj+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            // absorber at x=Nx grid boundary
            // Ex: odd-even-even
            EB_WAVE_ref(Nx-2-offset+1,jj  ,kk  ) = E_Xdir_OLD_ref(4+1,jj  ,kk  )
                                                 + cnst * ( EB_WAVE_ref(Nx-4-offset+1,jj  ,kk  ) 
                                                        -E_Xdir_OLD_ref(6+1          ,jj  ,kk  ) );
            // Ey: even-odd-even
            EB_WAVE_ref(Nx-2-offset  ,jj+1,kk  ) = E_Xdir_OLD_ref(4  ,jj+1,kk  )
                                                 + cnst * ( EB_WAVE_ref(Nx-4-offset  ,jj+1,kk  ) 
                                                        -E_Xdir_OLD_ref(6            ,jj+1,kk  ) );
            // Ez: even-even-odd
            EB_WAVE_ref(Nx-2-offset  ,jj  ,kk+1) = E_Xdir_OLD_ref(4  ,jj  ,kk+1)
                                                 + cnst * ( EB_WAVE_ref(Nx-4-offset  ,jj  ,kk+1) 
                                                        -E_Xdir_OLD_ref(6            ,jj  ,kk+1) );
        }
    }

    // absorber into y-direction
    //printf("abs_Mur_1st_v2: y1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            // absorber at y=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE_ref(ii+1,offset+0  ,kk  ) = E_Ydir_OLD_ref(ii+1,2  ,kk  )
                                              + cnst * ( EB_WAVE_ref(ii+1,offset+2  ,kk  )
                                                     -E_Ydir_OLD_ref(ii+1,0         ,kk  ) );
            // Ey: even-odd-even
            EB_WAVE_ref(ii  ,offset+0+1,kk  ) = E_Ydir_OLD_ref(ii  ,2+1,kk  )
                                              + cnst * ( EB_WAVE_ref(ii  ,offset+2+1,kk  )
                                                     -E_Ydir_OLD_ref(ii  ,0+1       ,kk  ) );
            // Ez: even-even-odd
            EB_WAVE_ref(ii  ,offset+0  ,kk+1) = E_Ydir_OLD_ref(ii  ,2  ,kk+1)
                                              + cnst * ( EB_WAVE_ref(ii  ,offset+2  ,kk+1)
                                                     -E_Ydir_OLD_ref(ii  ,0         ,kk+1) );
        }
    }
    
   
    //printf("abs_Mur_1st_v2: y2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
            // absorber at y=Ny grid boundary
            // Ex: odd-even-even
            EB_WAVE_ref(ii+1,Ny-2-offset  ,kk  ) = E_Ydir_OLD_ref(ii+1,4  ,kk  )
                                                 + cnst * ( EB_WAVE_ref(ii+1,Ny-4-offset  ,kk  )
                                                        -E_Ydir_OLD_ref(ii+1,6            ,kk  ) );
            // Ey: even-odd-even
            EB_WAVE_ref(ii  ,Ny-2-offset+1,kk  ) = E_Ydir_OLD_ref(ii  ,4+1,kk  )
                                                 + cnst * ( EB_WAVE_ref(ii  ,Ny-4-offset+1,kk  )
                                                        -E_Ydir_OLD_ref(ii  ,6+1          ,kk  ) );
            // Ez: even-even-odd
            EB_WAVE_ref(ii  ,Ny-2-offset  ,kk+1) = E_Ydir_OLD_ref(ii  ,4  ,kk+1)
                                                 + cnst * ( EB_WAVE_ref(ii  ,Ny-4-offset  ,kk+1)
                                                        -E_Ydir_OLD_ref(ii  ,6            ,kk+1) );
        }
    }

    // absorber into z-direction
    //printf("abs_Mur_1st_v2: z1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            // absorber at z=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE_ref(ii+1,jj  ,offset+0) = E_Zdir_OLD_ref(ii+1,jj  ,2  )
                                            + cnst * ( EB_WAVE_ref(ii+1,jj  ,offset+2  )
                                                   -E_Zdir_OLD_ref(ii+1,jj  ,0         ) );
            // Ey: even-odd-even
            EB_WAVE_ref(ii  ,jj+1,offset+0) = E_Zdir_OLD_ref(ii  ,jj+1,2  )
                                            + cnst * ( EB_WAVE_ref(ii  ,jj+1,offset+2  )
                                                   -E_Zdir_OLD_ref(ii  ,jj+1,0         ) );
            // Ez: even-even-odd
            EB_WAVE_ref(ii  ,jj  ,offset+0+1) = E_Zdir_OLD_ref(ii  ,jj  ,2+1)
                                              + cnst * ( EB_WAVE_ref(ii  ,jj  ,offset+2+1)
                                                     -E_Zdir_OLD_ref(ii  ,jj  ,0+1       ) );
        }
    }
    
    //printf("abs_Mur_1st_v2: z2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            // absorber at z=Nz grid boundary
            // Ex: odd-even-even
            EB_WAVE_ref(ii+1,jj  ,Nz_ref-2-offset  ) = E_Zdir_OLD_ref(ii+1,jj  ,4  )
                                                     + cnst * ( EB_WAVE_ref(ii+1,jj  ,Nz-4-offset  )
                                                            -E_Zdir_OLD_ref(ii+1,jj  ,6            ) );
            // Ey: even-odd-even
            EB_WAVE_ref(ii  ,jj+1,Nz_ref-2-offset  ) = E_Zdir_OLD_ref(ii  ,jj+1,4  )
                                                     + cnst * ( EB_WAVE_ref(ii  ,jj+1,Nz-4-offset  )
                                                            -E_Zdir_OLD_ref(ii  ,jj+1,6            ) );
            // Ez: even-even-odd
            EB_WAVE_ref(ii  ,jj  ,Nz_ref-2-offset+1) = E_Zdir_OLD_ref(ii  ,jj  ,4+1)
                                                     + cnst * ( EB_WAVE_ref(ii  ,jj  ,Nz-4-offset+1)
                                                            -E_Zdir_OLD_ref(ii  ,jj  ,6+1          ) );
        }
    }

    return EXIT_SUCCESS;

} //}}}

int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               systemGrid *G ) {
    //{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk;

    double
        aux, ny;

    ny  = 1e-4;
    
#pragma omp parallel for default(shared) private(ii,jj,kk,aux)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                aux = ny*(   EB_WAVE(ii  ,jj  ,kk+1+2) 
                          +  EB_WAVE(ii  ,jj  ,kk+1-2)
                          -2*EB_WAVE(ii  ,jj  ,kk+1  )
                        );
                EB_WAVE(ii  ,jj  ,kk+1) += aux;
            }
        }
    }

    return EXIT_SUCCESS;
}//}}}

/*UPML functions*/
double sigma(int pml_size, double nn, int m, double ds){

    double sig, sig_max, R_0;

    R_0 = pow(10,-6);
    sig_max = -(m+1)*log( R_0 )/(2*pml_size*ds);

    sig = pow( (nn) /(pml_size), m) * sig_max;

    return sig;  
}

void init_UPML_parameters(   gridConfiguration *gridCfg, boundaryGrid *boundaryG){

    int ii, jj, kk, count;
    double sig, kx, ky, kz;
    
    count = d_absorb;
    kx = 1;
    ky = 1;
    kz = 1;

    for ( ii=1 ; ii < (Nx/2)-1 ; ii+=1 ) {
        if(ii <= d_absorb + 1){

            sig = sigma(d_absorb, count, 4, dx);
            F1x(ii) = (2*kx) + (sig*dt);
            F2x(ii) = (2*kx) - (sig*dt);
            Cx(ii) = F2x(ii)/F1x(ii);

            count -= 1;
        }else if( ii >= (Nx/2) - d_absorb - 2){
            count += 1; 

            sig = sigma(d_absorb, count, 4, dx);
            F1x(ii) = (2*kx) + (sig*dt);
            F2x(ii) = (2*kx) - (sig*dt);
            Cx(ii) = F2x(ii)/F1x(ii);  
        }
    }

    for ( jj=1 ; jj < (Ny/2)-1 ; jj+=1 ) {
        if(jj <= d_absorb + 1){

            sig = sigma(d_absorb, count, 4, dx);
            F1y(jj) = (2*ky) + (sig*dt);
            F2y(jj) = (2*ky) - (sig*dt);
            Cy(jj) = F2y(jj)/F1y(jj);

            count -= 1;
        }else if( jj >= (Ny/2) - d_absorb - 2){
            count += 1; 

            sig = sigma(d_absorb, count, 4, dx);
            F1y(jj) = (2*ky) + (sig*dt);
            F2y(jj) = (2*ky) - (sig*dt);
            Cy(jj) = F2y(jj)/F1y(jj);  
        }
    }

    for ( kk=1 ; kk < (Nz/2)-1 ; kk+=1 ) {
        if(kk <= d_absorb + 1){

            sig = sigma(d_absorb, count, 4, dx);
            F1z(kk) = (2*kz) + (sig*dt);
            F2z(kk) = (2*kz) - (sig*dt);
            Cz(kk) = F2z(kk)/F1z(kk);

            count -= 1;
        }else if( kk >= (Nz/2) - d_absorb - 2){
            count += 1; 

            sig = sigma(d_absorb, count, 4, dx);
            F1z(kk) = (2*kz) + (sig*dt);
            F2z(kk) = (2*kz) - (sig*dt);
            Cz(kk) = F2z(kk)/F1z(kk); 
        }
    }

    for ( kk=1 ; kk < (Nz_ref/2)-1 ; kk+=1 ) {
        if(kk <= d_absorb + 1){

            sig = sigma(d_absorb, count, 4, dx);
            F1zr(kk) = (2*kz) + (sig*dt);
            F2zr(kk) = (2*kz) - (sig*dt);
            Czr(kk) = F2zr(kk)/F1zr(kk);

            count -= 1;
        }else if( kk >= (Nz_ref/2) - d_absorb - 2){
            count += 1; 

            sig = sigma(d_absorb, count, 4, dx);
            F1zr(kk) = (2*kz) + (sig*dt);
            F2zr(kk) = (2*kz) - (sig*dt);
            Czr(kk) = F2zr(kk)/F1zr(kk); 
        }
    }

}