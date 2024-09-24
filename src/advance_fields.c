#include "advance_fields.h"

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
    for (ii=2 ; ii < Nx-2 ; ii+=2) {
        for (jj=2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE(ii  ,jj+1,kk+1) += -1. * dt/dx * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                            );

                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE(ii+1,jj  ,kk+1) += -1. * dt/dx * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                            );

                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE(ii+1,jj+1,kk  ) += -1. * dt/dx * (
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
                EB_WAVE_ref(ii  ,jj+1,kk+1) += -1. * dt/dx * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                                );

                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE_ref(ii+1,jj  ,kk+1) += -1. * dt/dx * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                                );

                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE_ref(ii+1,jj+1,kk  ) += -1. * dt/dx * (
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

void advance_fields( gridConfiguration *gridCfg, 
                    systemGrid *G){
    
    //Advance wave-plasma current
    advance_J(      gridCfg, G );

    // advance wave magnetic field
    advance_B(      gridCfg, G );
    advance_B_ref(  gridCfg, G );

    // advance wave electric field
    advance_E(      gridCfg, G );
    advance_E_ref(  gridCfg, G );

}

