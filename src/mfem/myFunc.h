#ifndef MYFUNC_H_
#define MYFUNC_H_

#include "util.h"
//#include "eutectic.h"
#include "eutectic_rescaled.h"

// Boundary and initial physical attribute for mechanics

enum bndryType {dirichlet, neumann, robin, missed};

typedef struct {

    double theta ;
    double mu_s  ;
    double mu_f  ;
    double rho_f ;
    double rho_s ;
    double rho_r ;
    double gx    ;
    double gy    ;
    double invk0 ;
    double phi0  ;
    double U0    ;
    double phi_f_hat ;
    double V0    ;

    double l0    ;
    double u0    ;
    double p0    ;

    double L0    ;

    double l;

    // densities of two species
    double rho_1;
    double rho_2;

} PhysProperty;

// Extension class using phase package ==========================

class Phase {
    public:
        Phase() {};
        ~Phase() {delete pp; delete pPtr;};

        PhysProperty * pp;
        EUTECTIC::phase * pPtr;
};

// ==============================================================

double AssignPorosity(const vertex& point, PhysProperty * pp);

double AssignPorosity(double phi_f);

double AssignPorosity(const vertex& point, Phase * phase);

void AssignPhyProperties(PhysProperty * pp);
// ===================================================
// Define boundary condition
// ===================================================

// True solution
std::array<double, 3> trueSol(const vertex& point);

// Dirichlet value defined on the boundary
const vertex Dirichlet_val(const vertex& point);

// Return source term for darcy system as sum of velocity and pressure gradient  
const vertex darcyForce(const vertex& point, PhysProperty * pp);

const vertex darcyForce(const vertex& point, Phase * phase);

// Return source term for stokes system as sum of true solutions 
const vertex stokesForce(const vertex& point, PhysProperty * pp);

const vertex stokesForce(const vertex& point, Phase * phase);

// Return traction defined on the boundary
const vertex traction(const vertex& point, PhysProperty * pp);

const vertex traction(const vertex& point, Phase * phase);

// Boundary Condition
vertex bndryVs(const vertex& point, PhysProperty * pp);

vertex bndryVs(const vertex& point, Phase * phase);

vertex essenbndryVs(const vertex& point, PhysProperty * pp);
double naturbndryVs(const vertexSet& edge, PhysProperty * pp);

vertex bndryu(const vertex& point, PhysProperty * pp);

vertex bndryu(const vertex& point, Phase * phase);

double naturvalStokes(const vertex& point, PhysProperty * pp);

const bndryType bndryTypeMarker(const MeshInfo& mi, 
                                const indice& global,
                                const int& local);

const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global,
                                const int& local,
                                const std::vector<double>& parameter);

const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global);

const bndryType bndryTypeMarkerStokes(const MeshInfo& mi,
                                      const indice& global,
                                      const int& local);

const bndryType bndryTypeMarkerDarcy(const MeshInfo& mi,
                                      const indice& global,
                                      const int& edge);


#endif
