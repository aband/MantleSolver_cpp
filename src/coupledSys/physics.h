#ifndef PHYSICS_H_
#define PHYSICS_H_

// A struct object containing all physical attributes
typedef struct {

    double theta ;
    double mu_s  ;
    double mu_f  ;
    double rho_f ;
    double rho_s ;
    double gx    ;
    double gy    ;
    double invk0 ;
    double phi0  ;
    double U0    ;
    double phi_f_hat ;

    double l0    ;
    double u0    ;
    double p0    ;

    double L0    ;

    double l;

    // densities of two species
    double rho_1;
    double rho_2;

} PhysProperty;

#endif
