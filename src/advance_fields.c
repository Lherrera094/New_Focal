#include "advance_fields.h"

void advance_fields( gridConfiguration *gridCfg, 
                     systemGrid *G,
                     boundaryGrid *boundaryG){
    
    if(boundary_sel != 3){

        //Advance wave-plasma current
        advance_J(      gridCfg, G );

        // advance wave magnetic field
        advance_B(      gridCfg, G );
        advance_B_ref(  gridCfg, G );

        // advance wave electric field
        advance_E(      gridCfg, G );
        advance_E_ref(  gridCfg, G );

    }else if(boundary_sel == 3){

        //Advance wave-plasma current
        advance_J(      gridCfg, G );

        //Advance wave magnetic field
        advance_B_PML(  gridCfg, G, boundaryG );
        advance_Bref_PML(   gridCfg, G, boundaryG );

        //Advance wave electric field
        advance_E_PML(  gridCfg, G, boundaryG );
        advance_Eref_PML(  gridCfg, G, boundaryG );

    }
}

int advance_J( gridConfiguration *gridCfg, 
               systemGrid *G ) {
//{{{
    // This functions advances the current density J in time. J is calculated
    // from the fluid equation of motion of the electrons and reads
    // J_new = J_old + epsion_0*w_pe^2*E - w_ce*(Jx\hat(B)_0) - nu*J
    // Note that w_pe^2 --> n_e and w_ce --> B_0 with \hat(B) being the unit
    // vector pointing into the direction of B_0.
    // nu is a term corresponding to collisional damping that can be used to
    // respresent the effect of collisional and/or avoid some numerical 
    // instabilities that might arise at resonance like the upper-hybrid
    // resonance.

    size_t
        ii, jj, kk;

    // Jx: odd-even-even
    // Jy: even-odd-even
    // Jz: even-even-odd
    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                // Jx: odd-even-even
                J_B0(ii+1,jj  ,kk  ) += dt * (
                                        pow(2*M_PI,2) * n_e(ii/2,jj/2,kk/2) * EB_WAVE(ii+1,jj  ,kk  )
                                        - 2*M_PI * ( J_B0(ii  ,jj+1,kk  ) * J_B0(ii+1,jj+1,kk  )        // +Jy*B0z
                                                    -J_B0(ii  ,jj  ,kk+1) * J_B0(ii+1,jj  ,kk+1)        // -Jz*B0y
                                                    )
                                            );

                // Jy: even-odd-even
                J_B0(ii  ,jj+1,kk  ) += dt * (
                                        pow(2*M_PI,2) * n_e(ii/2,jj/2,kk/2) * EB_WAVE(ii  ,jj+1,kk  )
                                        - 2*M_PI * (-J_B0(ii+1,jj  ,kk  ) * J_B0(ii+1,jj+1,kk  )         // -Jx*B0z
                                                    +J_B0(ii  ,jj  ,kk+1) * J_B0(ii  ,jj+1,kk+1)         // +Jz*B0x
                                                    )
                                            );

                // Jz: even-even-odd
                J_B0(ii  ,jj  ,kk+1)  += dt * (
                                        pow(2*M_PI,2) * n_e(ii/2,jj/2,kk/2) * EB_WAVE(ii  ,jj  ,kk+1)
                                        - 2*M_PI * ( J_B0(ii+1,jj  ,kk  ) * J_B0(ii+1,jj  ,kk+1)        // +Jx*B0y
                                                    -J_B0(ii  ,jj+1,kk  ) * J_B0(ii  ,jj+1,kk+1)        // -Jy*B0x
                                                    )
                                            );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_B( gridConfiguration *gridCfg, 
               systemGrid *G ) {
//{{{
    // B_new = B_old - nabla x E

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=4 ; ii < Nx-2 ; ii+=2) {
        for (jj=4 ; jj < Ny-2 ; jj+=2) {
            for (kk=4 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE(ii  ,jj+1,kk+1) += -1. * (dt/dx) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                            );

                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE(ii+1,jj  ,kk+1) += -1. * (dt/dx) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                            );

                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE(ii+1,jj+1,kk  ) += -1. * (dt/dx) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                            );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_B_ref( gridConfiguration *gridCfg, 
                   systemGrid *G ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE_ref(ii  ,jj+1,kk+1) += -1. * (dt/dx) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                                );

                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE_ref(ii+1,jj  ,kk+1) += -1. * (dt/dx) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                                );

                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE_ref(ii+1,jj+1,kk  ) += -1. * (dt/dx) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                                );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E( gridConfiguration *gridCfg, 
               systemGrid *G ) {
//{{{
    // E_new = E_old + c^2*nablaxB - 1/epsilon_0*J

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE(ii+1,jj  ,kk  ) += dt/dx * (
                                        +EB_WAVE(ii+1,jj+1,kk  ) - EB_WAVE(ii+1,jj-1,kk  )
                                        -EB_WAVE(ii+1,jj  ,kk+1) + EB_WAVE(ii+1,jj  ,kk-1)
                                        ) - dt * J_B0(ii+1,jj  ,kk  );

                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE(ii  ,jj+1,kk  ) += dt/dx * (
                                        +EB_WAVE(ii  ,jj+1,kk+1) - EB_WAVE(ii  ,jj+1,kk-1)
                                        -EB_WAVE(ii+1,jj+1,kk  ) + EB_WAVE(ii-1,jj+1,kk  )
                                        ) - dt * J_B0(ii  ,jj+1,kk  );

                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE(ii  ,jj  ,kk+1) += dt/dx * (
                                        +EB_WAVE(ii+1,jj  ,kk+1) - EB_WAVE(ii-1,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+1) + EB_WAVE(ii  ,jj-1,kk+1)
                                        ) - dt * J_B0(ii  ,jj  ,kk+1);
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E_ref( gridConfiguration *gridCfg, 
                   systemGrid *G ) { 
//{{{
    // same as advance_E but for reference fields (directional coupler)
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE_ref(ii+1,jj  ,kk  ) += dt/dx * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE_ref(ii  ,jj+1,kk  ) += dt/dx * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE_ref(ii  ,jj  ,kk+1) += dt/dx * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


/*UPML Functions*/
int advance_B_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ) {
//{{{      

    size_t
        ii, jj, kk;
    /*double
        dxstore, dystore, dzstore;*/

/*Core grid advancing*/
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE(ii  ,jj+1,kk+1) += -1. * (dt/dx) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE(ii+1,jj  ,kk+1) += -1. * (dt/dx) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE(ii+1,jj+1,kk  ) += -1. * (dt/dx) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                // -dBx/dt = dEz/dy - dEy/dz
                /*dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore 
                                        );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore);

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore);*/
                
            }
        }
    }

    UPML_B_faces(  gridCfg, G, boundaryG );
    UPML_B_corners(gridCfg, G, boundaryG );
    UPML_B_edges(  gridCfg, G, boundaryG );

    return EXIT_SUCCESS;
}//}}}

int advance_Bref_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ) {
//{{{      

    size_t
        ii, jj, kk;


/*Core grid advancing*/
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE_ref(ii  ,jj+1,kk+1) += -1. * (dt/dx) * (
                                        +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                        -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                        );

                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE_ref(ii+1,jj  ,kk+1) += -1. * (dt/dx) * (
                                        +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                        -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                        );

                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE_ref(ii+1,jj+1,kk  ) += -1. * (dt/dx) * (
                                        +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                        -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                        );
                
            }
        }
    }

    UPML_Bref_faces(  gridCfg, G, boundaryG );
    //UPML_Bref_corners(gridCfg, G, boundaryG);
    //UPML_Bref_edges(  gridCfg, G, boundaryG );

    return EXIT_SUCCESS;
}//}}}

int advance_E_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ){
//{{{
    size_t
        ii, jj, kk;
    /*double
        dxstore, dystore, dzstore;*/

/*Core grid advancing*/
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE(ii+1,jj  ,kk  ) += dt/dx * (
                                        +EB_WAVE(ii+1,jj+1,kk  ) - EB_WAVE(ii+1,jj-1,kk  )
                                        -EB_WAVE(ii+1,jj  ,kk+1) + EB_WAVE(ii+1,jj  ,kk-1)
                                        ) - dt * J_B0(ii+1,jj  ,kk  );

                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE(ii  ,jj+1,kk  ) += dt/dx * (
                                        +EB_WAVE(ii  ,jj+1,kk+1) - EB_WAVE(ii  ,jj+1,kk-1)
                                        -EB_WAVE(ii+1,jj+1,kk  ) + EB_WAVE(ii-1,jj+1,kk  )
                                        ) - dt * J_B0(ii  ,jj+1,kk  );

                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE(ii  ,jj  ,kk+1) += dt/dx * (
                                        +EB_WAVE(ii+1,jj  ,kk+1) - EB_WAVE(ii-1,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+1) + EB_WAVE(ii  ,jj-1,kk+1)
                                        ) - dt * J_B0(ii  ,jj  ,kk+1);

                // dEx/dt = (dBz/dy - dBy/dz)
                /*dxstore = DH_WAVE(ii+1,jj  ,kk  );
                DH_WAVE(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii+1,jj+1,kk  ) - EB_WAVE(ii+1,jj-1,kk  )
                                        -EB_WAVE(ii+1,jj  ,kk+1) + EB_WAVE(ii+1,jj  ,kk-1)
                                        );

                EB_WAVE(ii+1,jj  ,kk  ) = Cz(kk/2)*EB_WAVE(ii+1,jj  ,kk  ) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE(ii  ,jj+1,kk  );
                DH_WAVE(ii  ,jj+1,kk  ) = Cz(kk/2)*DH_WAVE(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii  ,jj+1,kk+1) - EB_WAVE(ii  ,jj+1,kk-1)
                                        -EB_WAVE(ii+1,jj+1,kk  ) + EB_WAVE(ii-1,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE(ii  ,jj  ,kk+1);
                DH_WAVE(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+1) - EB_WAVE(ii-1,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+1) + EB_WAVE(ii  ,jj-1,kk+1)
                                        );

                EB_WAVE(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );*/
            }
        }
    }

    UPML_E_faces(  gridCfg, G, boundaryG );
    UPML_E_corners(gridCfg, G, boundaryG );
    UPML_E_edges(  gridCfg, G, boundaryG );

    return EXIT_SUCCESS;
}


int advance_Eref_PML(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ){
//{{{
    size_t
        ii, jj, kk;

/*Core grid advancing*/
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE_ref(ii+1,jj  ,kk  ) += dt/dx * (
                                        +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                        -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                        );

                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE_ref(ii  ,jj+1,kk  ) += dt/dx * (
                                        +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                        -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                        );

                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE_ref(ii  ,jj  ,kk+1) += dt/dx * (
                                        +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                        -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                        );
            }
        }
    }

    //UPML_Eref_faces(  gridCfg, G, boundaryG );
    //UPML_Eref_corners(gridCfg, G, boundaryG );
    //UPML_Eref_edges(  gridCfg, G, boundaryG );

    return EXIT_SUCCESS;
}
