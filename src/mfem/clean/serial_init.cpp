#include "serial_solver.h"

int DarcyStokes::init(const MeshInfo& mi, PhysProperty * pp, const std::vector<double>& param){

    basis_ = basis();
    hdiv_  = Hdivmixed();
    br_    = BRMixed();

    br_.ComputeTotalDOF(mi);
    hdiv_.ComputeTotalDOF(mi);

    ComputeEssenBndryAll(mi,pp,param);

    return 1;
}

int DarcyStokes::printBndryAll(){

    // All the assigned values for essential Stokes boundary values
    for (auto& it: bndryStokesEssenAll){

        cout << "Global dof : " << it->first << it-><< ; 


    }

    // All the assigned values for essential Darcy boundary values
    for (auto& it: bndryDarcyEssenAll){


    }

    return 1;
}
