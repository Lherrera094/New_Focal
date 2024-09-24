#ifndef MACROS_BOUNDARY_H
#define MACROS_BOUNDARY_H

#define ABSORBER_DAMPING(eco,damp) (1.-eco*damp*damp)

/*Macros for ABC*/
#define ecoBG(boundaryG)            boundaryG->eco

#define eco                         ecoBG(boundaryG)

/*Macros MUR abosbing boundary*/


/*Macros for PML*/
      

#endif
