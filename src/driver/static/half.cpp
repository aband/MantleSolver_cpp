#include "myFunc.h"

void AssignPhyProperties(PhysProperty * pp){

    pp->theta = 0.0;
    pp->mu_s  = 1e15;
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
}


double AssignPorosity(const vertex& point, PhysProperty * pp){

    if (abs(point[1]) < 120*1000/pp->l0 && abs(point[0]) < abs(point[1]) + pp->l){

        double value = 0.05*pow((120*1000/pp->l0 - abs(point[1]))/(120*1000/pp->l0),2) * 
                               (1-abs(point[0])/(abs(point[1])+pp->l));

        //return value;
        return 0.0;
    } else {
        return 0.0;
    }

}

double AssignPorosity(double phi_f){
    return 1-phi_f;
}

// =============================================================================

// Boundary values
// Constant upwelling velocity ascending model
vertex bndryVs(const vertex& point, PhysProperty * pp){

/*
    // Stokes
    double V0 = pp->V0 / pp->u0 * -1;

    vertex work = {0.0,0.0};

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

    return work; 
*/
    double V0 = pp->V0 / pp->u0 * -1;
/*
    if (point[1] < -0.1){
        if (point[0]<0.05){
            //return {0.0, -V0 * sin(point[0]/0.05*3.14159265358979323846/2.0)};
				return {0.0, -V0};
        }else{
            return {0.0, -V0}; 
        }
    } else {
        return {0.0,0.0};
    }
*/

    //double velhead = -1*V0;
    double velhead = 1;

    vertex work {0.0,0.0};

    if (point[1] > -0.0005){
        if (point[0]<0.05){
            work =  {velhead * sin(point[0]/0.05*3.14159265358979323846/2.0), 0.0};
				//return {1, 0.0};
        }else{
            work = {velhead, 0.0}; 
        }
    } else if (point[1] < -0.495){
 
        work = {0.0,velhead};

    } else {
        work = {0.0,0.0};
    }

    return work; 


//    return {0.0, -V0};

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
    double phi_f = AssignPorosity(point, pp);

    double rho_r = pp->rho_f*phi_f + pp->rho_s*(1-phi_f);

    //double coef1 = (1-pp->phi0) * pow(pp->phi0,2+2*pp->theta);
    //double coef2 = 4*pp->mu_s*pp->U0/(3.14159265358979323846*(x*x+z*z)*pp->x0*pp->x0) /rho_r /pp->gy;
    double coef1 = (1-phi_f)*pow(phi_f,2+2*pp->theta); 

    double coef2 = 4*pp->U0/pp->u0/(3.14159265358979323846*(x*x+z*z)*(x*x+z*z));

    work[0] = coef1*coef2*2*x*z;
    work[1] = coef1*coef2*(z*z-x*x);

    work[0] += coef1 * 0;
    work[1] += coef1 * 1;

    return {0.0, 0.0};

    //return work;
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


    if (point[1]<-0.49){
        return {0.0,0.0};
    } else {
        return {0.0,-1.0}; 
    }


//    return {0.0, 0.0};
}

const vertex traction(const vertex& point, PhysProperty * pp){
    // return traction defined on the boundary
    // zero traction situation

    if (point[0] > 0.1){

        return {1.0,0.0}; // free stress

    } else {
        return {0.0,0.0};
    }

}

// Mark boundary type for stokes
/*
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

    // Top edge all normal component are set to be zero
    // No flow out of the domain on top
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        //it = top_normal.find(local);
        it = top_tang.find(local);
        if (it != top_tang.end()){
            type = neumann;
        }
		  //type = dirichlet;
    } 

    // left side symmetrical condition
    // no normal flux

    // right side free outflow
    if (global[0] == 0) {
        it = left_tang.find(local);
        if (it != left_tang.end()){
            type = neumann;
        }
		  //type = dirichlet;
    } 
 
    if (global[0] == mi.MPIglobalCellSize[0]-1){
//        it = right_tang.find(local);
//        if (it != right_tang.end()){
//            type = neumann;
//        }
        type = neumann;
        //type = dirichlet;
    }

    if (global[1] == 0 ){
        type = dirichlet;
    }

    return type;
}
*/

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

    // All dirichlet at top
    if (global[1] == mi.MPIglobalCellSize[1]-1){
        //it = top_normal.find(local);
        //it = top_tang.find(local);
        //if (it != top_tang.end()){
        //    type = neumann;
        //}
        type = dirichlet;
    } 

    // left side symmetrical condition
    // no normal flux

    // right side free outflow
    if (global[0] == 0) {
        it = left_tang.find(local);
        if (it != left_tang.end()){
            type = neumann;
        }
		  //type = dirichlet;
    } 
 
    if (global[0] == mi.MPIglobalCellSize[0]-1){
//        it = right_tang.find(local);
//        if (it != right_tang.end()){
//            type = neumann;
//        }
        type = neumann;
        //type = dirichlet;
    }

    if (global[1] == 0 ){

        if (global[0] < mi.MPIglobalCellSize[0]/6*5){
            type = neumann;  
        }else {
            it = bottom_normal.find(local);
            if (it != bottom_normal.end()) {
                type = dirichlet;
            } else {
                type = neumann;
            }
        }
    }

    return type;
}


// Mark boundary type for darcy
const bndryType bndryTypeMarkerDarcy(const MeshInfo& mi,
                                     const indice& global,
                                     const int& edge){

    return dirichlet;
}
