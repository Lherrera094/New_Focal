#include "UPML_module.h"

/*Magnetic field UPML*/
void UPML_B_faces(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ){
    
    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_B_corners(gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb+2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb+2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb+2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb+2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_B_edges(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Edge x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );    
            }
        }
    }

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );    
            }
        }
    }

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );    
            }
        }
    }

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE(ii  ,jj+1,kk+1);
                DH_WAVE(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                        +EB_WAVE(ii  ,jj+2,kk+1) - EB_WAVE(ii  ,jj  ,kk+1)
                                        -EB_WAVE(ii  ,jj+1,kk+2) + EB_WAVE(ii  ,jj+1,kk  )
                                        );

                EB_WAVE(ii  ,jj+1,kk+1) = Cz(kk/2)*EB_WAVE(ii  ,jj+1,kk+1) + ( 1/F2z(kk/2) )*(
                                        + F2x(ii/2)*DH_WAVE(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE(ii+1,jj  ,kk+1);
                DH_WAVE(ii+1,jj  ,kk+1) = Cz(kk/2)*DH_WAVE(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2z(kk/2) ) * (
                                        +EB_WAVE(ii+1,jj  ,kk+2) - EB_WAVE(ii+1,jj  ,kk  )
                                        -EB_WAVE(ii+2,jj  ,kk+1) + EB_WAVE(ii  ,jj  ,kk+1)
                                        );

                EB_WAVE(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                        + F2y(jj/2)*DH_WAVE(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE(ii+1,jj+1,kk  );
                DH_WAVE(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                        +EB_WAVE(ii+2,jj+1,kk  ) - EB_WAVE(ii  ,jj+1,kk  )
                                        -EB_WAVE(ii+1,jj+2,kk  ) + EB_WAVE(ii+1,jj  ,kk  )
                                        );

                EB_WAVE(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                        + F2z(kk/2)*DH_WAVE(ii+1,jj+1,kk  ) - F1z(kk/2)*dzstore );    
            }
        }
    }

}

void UPML_Bref_faces(   gridConfiguration *gridCfg, 
                        systemGrid *G,
                        boundaryGrid *boundaryG ){
    
    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Bref_corners(gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb+2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb+2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb+2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb+2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Bref_edges(   gridConfiguration *gridCfg, 
                        systemGrid *G,
                        boundaryGrid *boundaryG ){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Edge x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );   
            }
        }
    }

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );   
            }
        }
    }

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );    
            }
        }
    }

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref(ii  ,jj+1,kk+1);
                DH_WAVE_ref(ii  ,jj+1,kk+1) = Cy(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - 1.*( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+2,kk+1) - EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+2) + EB_WAVE_ref(ii  ,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk+1) = Czr(kk/2)*EB_WAVE_ref(ii  ,jj+1,kk+1) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii  ,jj+1,kk+1) - F1x(ii/2)*dxstore );

                // -dBy/dt = dEx/dz - dEz/dx
                dystore = DH_WAVE_ref(ii+1,jj  ,kk+1);
                DH_WAVE_ref(ii+1,jj  ,kk+1) = Czr(kk/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - 1.*( 2*dt/dx/F2zr(kk/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+2) - EB_WAVE_ref(ii+1,jj  ,kk  )
                                            -EB_WAVE_ref(ii+2,jj  ,kk+1) + EB_WAVE_ref(ii  ,jj  ,kk+1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk+1) = Cx(ii/2)*EB_WAVE_ref(ii+1,jj  ,kk+1) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk+1) - F1y(jj/2)*dystore );

                // -dBz/dt = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref(ii+1,jj+1,kk  );
                DH_WAVE_ref(ii+1,jj+1,kk  ) = Cx(ii/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - 1.*( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+2,jj+1,kk  ) - EB_WAVE_ref(ii  ,jj+1,kk  )
                                            -EB_WAVE_ref(ii+1,jj+2,kk  ) + EB_WAVE_ref(ii+1,jj  ,kk  )
                                            );

                EB_WAVE_ref(ii+1,jj+1,kk  ) = Cy(jj/2)*EB_WAVE_ref(ii+1,jj+1,kk  ) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii+1,jj+1,kk  ) - F1zr(kk/2)*dzstore );    
            }
        }
    }

}


/*Electric field UPML*/
void UPML_E_faces(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = Nx - d_absorb - 2 ; ii < Nx - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_E_corners(gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb+2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_E_edges(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ){

    // DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );   
            }
        }
    }

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz - d_absorb - 2 ; kk < Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE(ii+1,jj  ,kk  );
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
                                        + F2z(kk/2)*DH_WAVE(ii  ,jj  ,kk+1) - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Eref_faces(  gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG ){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = Nx - d_absorb - 2 ; ii < Nx - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj=Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Eref_corners(gridConfiguration *gridCfg, 
                    systemGrid *G,
                    boundaryGrid *boundaryG){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb+2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb+2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb+2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb+2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Eref_edges(   gridConfiguration *gridCfg, 
                        systemGrid *G,
                        boundaryGrid *boundaryG ){

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Edge x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = d_absorb + 2 ; kk < Nz_ref - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb + 2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = Nx - d_absorb - 2 ; ii < Nx-2 ; ii+=2) {                
        for (jj = d_absorb + 2 ; jj < Ny - d_absorb - 2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );   
            }
        }
    }

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb + 2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );  
            }
        }
    }

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb + 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );    
            }
        }
    }

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb + 2 ; ii < Nx - d_absorb - 2 ; ii+=2) {                
        for (jj = Ny - d_absorb - 2 ; jj < Ny-2 ; jj+=2) {
            for (kk = Nz_ref - d_absorb - 2 ; kk < Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref(ii+1,jj  ,kk  );
                DH_WAVE_ref(ii+1,jj  ,kk  ) = Cy(jj/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) + ( 2*dt/dx/F2y(jj/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj+1,kk  ) - EB_WAVE_ref(ii+1,jj-1,kk  )
                                            -EB_WAVE_ref(ii+1,jj  ,kk+1) + EB_WAVE_ref(ii+1,jj  ,kk-1)
                                            );

                EB_WAVE_ref(ii+1,jj  ,kk  ) = Czr(kk/2)*EB_WAVE_ref(ii+1,jj  ,kk  ) + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref(ii+1,jj  ,kk  ) - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref(ii  ,jj+1,kk  );
                DH_WAVE_ref(ii  ,jj+1,kk  ) = Czr(kk/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) + ( 2*dt/dx/F2z(kk/2) ) * (
                                            +EB_WAVE_ref(ii  ,jj+1,kk+1) - EB_WAVE_ref(ii  ,jj+1,kk-1)
                                            -EB_WAVE_ref(ii+1,jj+1,kk  ) + EB_WAVE_ref(ii-1,jj+1,kk  )
                                            );

                EB_WAVE_ref(ii  ,jj+1,kk  ) = Cx(ii/2)*EB_WAVE_ref(ii  ,jj+1,kk  ) + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref(ii  ,jj+1,kk  ) - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref(ii  ,jj  ,kk+1);
                DH_WAVE_ref(ii  ,jj  ,kk+1) = Cx(ii/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) + ( 2*dt/dx/F2x(ii/2) ) * (
                                            +EB_WAVE_ref(ii+1,jj  ,kk+1) - EB_WAVE_ref(ii-1,jj  ,kk+1)
                                            -EB_WAVE_ref(ii  ,jj+1,kk+1) + EB_WAVE_ref(ii  ,jj-1,kk+1)
                                            );

                EB_WAVE_ref(ii  ,jj  ,kk+1) = Cy(jj/2)*EB_WAVE_ref(ii  ,jj  ,kk+1) + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref(ii  ,jj  ,kk+1) - F1zr(kk/2)*dzstore );    
            }
        }
    }

}
