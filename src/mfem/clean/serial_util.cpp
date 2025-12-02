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
	 cout << "Essential boundary condition for Stokes." << endl;
    for (auto& it: bndryStokesEssenAll){

        printf("g dof : %d , l dof : %d, g cell : (%d, %d), essen val : %e \n", 
               it.first, it.second.localDOF, it.second.globalElem[0], 
               it.second.globalElem[1], it.second.essenval);

    }

    cout << endl;

    // All the assigned values for essential Darcy boundary values
    cout << "Essential boundary condition for Darcy." << endl;
    for (auto& it: bndryDarcyEssenAll){

        cout << "Global dof : " << it.first << 
                " local dof : " << it.second.localDOF << endl;

    }

    return 1;
}
