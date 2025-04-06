#include "myFunc.h"

void AssignPhyProperties(PhysProperty * pp){

    pp->theta = 0.0;
    pp->mu_s  = 1e19;
    pp->mu_f  = 1.0;
    pp->rho_f = 2800;
    pp->rho_s = 3300;
    pp->gx    = 0.0;
    pp->gy    = -10.0;
    pp->invk0 = 1.0/(1e-8);
    pp->phi0  = 0.4;
    pp->U0    = 1e-9;
    pp->L0    = 160*1000;
    pp->V0    = 3.2/100/(365*24*60*60); //3.2 (cm/y)

    // Non dimensionalization parameters

    pp->rho_r = pp->rho_s - pp->rho_f;

    pp->l0    = pow(pp->mu_s/pp->invk0/pp->mu_f,0.5);
    pp->p0    = pp->gy*pp->l0*pp->rho_r;
    pp->u0    = pp->gy*pp->rho_r/pp->mu_f/pp->invk0;

    pp->l = 20/pp->l0;
}

double AssignPorosity(const vertex& point, PhysProperty * pp){

    // used to identify incorrect porosity
    if (point[1] > -0.2){
        return 0.04;
    }else {
        return 2.3e-40;
    }
}

double AssignPorosity(double phi_f){
    return 1-phi_f;
}

// ===================================================
// Below not needed if no true solution posted

// =============================================================

// Boundary values
// Constant upwelling velocity ascending model
vertex bndryVs(const vertex& point, PhysProperty * pp){

    // Stokes
    double V0 = pp->V0 / pp->u0 * -1;

    return {0.0,V0};
    //return {0.0, 1.0};
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    // Darcy

    return {0.0,0.0};
}

// ==============================================================
const vertex darcyForce(const vertex& point, PhysProperty * pp){

    // Constant porosity.
    return {0.0,0.0};
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    // Returns nondimensionalized gravity.
    // Attention!!! It should not be scaled by porosity
	 // porosity scale will be added in another function
double V0 = pp->V0 / pp->u0;	
    //return {0.0, -1.0/V0};
    return {0.0,-1.0};
}

const vertex traction(const vertex& point, PhysProperty * pp){
    // return traction defined on the boundary
    // zero traction situation

    return {0.0,0.0}; // free stress
}

// =========================================================================
const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global,
                                const int& local,
                                const std::vector<double>& parameter){

    bndryType type = missed;

    // normal dof 
    std::set<int> top_normal {11,7,6};
    std::set<int> left_normal {0,3,8};
    std::set<int> right_normal {1,2,10};
    std::set<int> bottom_normal {4,5,9};

    // tangent dof
    std::set<int> top_tang {3,2};
    std::set<int> left_tang {7,4};
    std::set<int> right_tang {5,6};
    std::set<int> bottom_tang {0,1};

    std::set<int>::iterator it;

    type = dirichlet;

    // Top edge all normal component are set free 
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        it = top_normal.find(local);
        if (it != top_normal.end()){
            type = neumann;
        }
    } 


    if (global[0] == 0) {
        it = left_tang.find(local);
        if (it != left_tang.end()){
            type = neumann;
        }
    } 
 
    if (global[0] == mi.MPIglobalCellSize[0]-1){
        it = right_tang.find(local);
        if (it != right_tang.end()){
            type = neumann;
        }
    }


    if (global[1] == 0 ){
        type = dirichlet;
    }


    return type;
    //return dirichlet;
}

inline int getdofset(std::set<int>& directionSet, 
                     const MeshInfo& mi,
                     const std::string& direction){

    return 0;
}

// A more direct way of marking boundary type
// Marking boundary type with global dof index
const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const int& global){

    bndryType type = missed;

    std::set<int> bottom;
    std::set<int> top;
    std::set<int> left;
    std::set<int> right;

    return type;
}

bool exit(double range, int M, int i){

    if (i < M/2 + range && i > M/2 - range){
        return true;
    } else {
        return false;
    }

}

const bndryType bndryTypeMarkerDarcy(const MeshInfo& mi,
                                     const indice& global,
                                     const int& edge){

    bndryType type = missed;

    //if (global[1] == mi.MPIglobalCellSize[1]-1 && exit(1.05, mi.MPIglobalCellSize[0], global[0]) ){
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        if (edge == 3){
            type = neumann; 
        }else {
            type = dirichlet;
        }
    } else if (global[1] == 0){
        if (edge == 3){
            type = neumann;
        } else {
            type = dirichlet;
        }
    } else {
        type = dirichlet;
    }

    return type;
}
