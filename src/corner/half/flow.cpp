#include "myFunc.h"

void AssignPhyProperties(PhysProperty * pp){

    pp->theta = 0.0;
    pp->mu_s  = 1e19;
    pp->mu_f  = 1.0;
    pp->rho_f = 2800;
    pp->rho_s = 3300;
    pp->gx    = 0.0;
    pp->gy    = 10.0;
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

    cout << pp->l0 << endl;

}

// Essential boundary Boundary values =================================
// Constant upwelling velocity ascending model
vertex bndryVs(const vertex& point, PhysProperty * pp){

    double V0 = pp->V0 / pp->u0;

    vertex work = {0.0,0.0};

/*
    // Test case 3:
    // Corner Flow
    double x, z;

    if (point[0] < 0.0) {
        x = point[0] - pp->l;
    }else{
        x = point[0] + pp->l;
    }

    z = point[1];

    double coef = 2*pp->U0/(3.14159265358979323846*(x*x+z*z))/pp->u0;
    //coef = 2/(3.14159265358979323846*(x*x+z*z));

    work =  {atan(x/(-1*z))*(x*x+z*z) + x*z, z*z};

    work *= coef;

    if (z == 0){
        if (point[0] < 0){
            work[0] = -1*pp->U0/pp->u0;
        } else {
            work[0] = pp->U0/pp->u0;
        }
    }
*/


    double cut = 0.1;

    if (point[1] > -0.00005){

       if (point[0]<cut){
            work[0] = V0 * sin(point[0]/cut*3.14159265358979323846/2.0);
        } else {
            work[0] = V0;
        }

    } else if (point[1] < -0.00005){

       work[1] = V0;

    }


    return work; 
}

vertex bndryu(const vertex& point, PhysProperty * pp){

    // Darcy
    vertex work {0.0,0.0};
 
    double x,z;

    if (point[0] < 0.0){
        x = point[0] - pp->l;
    }else {
        x = point[0] + pp->l;
    }

    z = point[1];

    // Point wise porosity
    //double phi_f = AssignPorosity(point, pp);
    double phi_f = 0.0;

    double rho_r = pp->rho_f*phi_f + pp->rho_s*(1-phi_f);

    //double coef1 = (1-pp->phi0) * pow(pp->phi0,2+2*pp->theta);
    //double coef2 = 4*pp->mu_s*pp->U0/(3.14159265358979323846*(x*x+z*z)*pp->x0*pp->x0) /rho_r /pp->gy;
    double coef1 = (1-phi_f)*pow(phi_f,2+2*pp->theta); 

    double coef2 = 4*pp->U0/pp->u0/(3.14159265358979323846*(x*x+z*z)*(x*x+z*z));

    work[0] = coef1*coef2*2*x*z;
    work[1] = coef1*coef2*(z*z-x*x);

    work[0] += coef1 * 0;
    work[1] += coef1 * 1;

    return {0.0,0.0};
}

// ====================================================================
const vertex darcyForce(const vertex& point, PhysProperty * pp){

    // Constant porosity.
    return {0.0,0.0};
}

const vertex stokesForce(const vertex& point, PhysProperty * pp){

    // Returns nondimensionalized gravity.
    // Attention!!! It should not be scaled by porosity
	 // porosity scale will be added in another function
double V0 = pp->V0 / pp->u0;	
    return {0.0, -1.0};
    //return {0.0,0.0};
}

double naturvalStokes(const vertex& point, PhysProperty * pp){

    return abs(point[1]);
}

// ====================================================================
const bndryType bndryTypeMarker(const MeshInfo& mi,
                                const indice& global,
                                const int& local,
                                const std::vector<double>& parameter){

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

    bndryType type = dirichlet;

    // Symmetrical boundary condition 
    if (global[0] == 0) {

        if (global[1]<mi.MPIglobalCellSize[1]-5){

        it = left_tang.find(local);
        if (it != left_tang.end()){
            type = neumann;
        }

        } 
		  //else {
        
		  //type = dirichlet;
        
		  //}
		  
		  //type = dirichlet;
    } 

    // Right side free outflow
    if (global[0] == mi.MPIglobalCellSize[0]-1){

        //if (global[1] > mi.MPIglobalCellSize[1]-5){
        //    type = dirichlet;
        //} else {
        //   type = neumann;
        //}

        type = neumann;
    }
 
    // Bottom edge fixed inflow
    if (global[1] == 0 ){
        //it = bottom_tang.find(local);
        //if (it != top_tang.end()){
        //    type = neumann;
        //}

        type = dirichlet;
    }

    // Top edge dirichlet
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        //it = top_normal.find(local);
        //if (it != top_normal.end()){
        //    type = neumann;
        //}
        type = dirichlet;
    } 

    // specificaly treat top right corner dof
    if (global[0] == mi.MPIglobalCellSize[0]-1 && 
        global[1] == mi.MPIglobalCellSize[1]-1){

        if (local == 2){
            type = dirichlet;
        } else if (local == 6){
            type = neumann;
        } else if (local == 10){
            type = neumann;
        }

    }

    return type;
}

const bndryType bndryTypeMarkerDarcy(const MeshInfo& mi,
                                     const indice& global,
                                     const int& edge){

    return dirichlet;
}
