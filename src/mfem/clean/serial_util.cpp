#include "serial_solver.h"

int DarcyStokes::init(const MeshInfo& mi, PhysProperty * pp, const std::vector<double>& param){

    basis_ = basis();
    hdiv_  = Hdivmixed();
    br_    = BRMixed();

    br_.ComputeTotalDOF(mi);
    hdiv_.ComputeTotalDOF(mi);

    // Compute boundary conditions
    ComputeEssenBndryAll(mi,pp,param);

    // Mark different types of boundary dofs
    refArrayStokesEssen_ = new int[br_.getDOF()];
    CreateRefMap(br_, mi, refArrayStokesEssen_, refArrayStokesNatur_, bndryDOFStokesEssen_, bndryDOFStokesNatur_, {0.0});

    refArrayDarcyEssen_ = new int[hdiv_.getDOF()];
    CreateRefMap(hdiv_, mi, refArrayDarcyEssen_, refArrayDarcyNatur_, bndryDOFDarcyEssen_, bndryDOFDarcyNatur_, {0.0});

    return 1;
}

int DarcyStokes::printBndryAll(){

    // Make sure all the boundary information are computed correctly
    printf("essential count Stokes: %d, Darcy: %d, natural count Stokes: %d, Darcy: %d \n", bndryDOFStokesEssen_, bndryDOFDarcyEssen_, bndryDOFStokesNatur_, bndryDOFDarcyNatur_);

    // All the assigned values for essential Stokes boundary values
	 cout << "Essential boundary condition for Stokes." << endl;
    for (auto& it: bndryStokesAll){

        printf("g dof : %d , l dof : %d, g cell : (%d, %d), essen val : %e , natur val : %e \n", 
               it.first, it.second.localDOF, it.second.globalElem[0], 
               it.second.globalElem[1], it.second.essenval, it.second.naturval);

    }

    cout << endl;

    // All the assigned values for essential Darcy boundary values
    cout << "Essential boundary condition for Darcy." << endl;
    for (auto& it: bndryDarcyAll){

        printf("g dof : %d , l dof : %d, g cell : (%d, %d), essen val : %e , natur val : %e \n", 
               it.first, it.second.localDOF, it.second.globalElem[0], 
               it.second.globalElem[1], it.second.essenval, it.second.naturval);

    }

    return 1;
}

static int printRedSys(ReducedSys& redsys){

    // Print a reduced system
    cout << "Matrix M : " << endl;
    MatView(redsys.M, PETSC_VIEWER_STDOUT_WORLD);

    cout << "Essen boundary : " << endl;    
    VecView(redsys.g, PETSC_VIEWER_STDOUT_WORLD);

    return 1;
}

int DarcyStokes::showMatrix(){

    printRedSys(reducedDarcy_); 

    return 1;
}
