#include "serial_solver.h"

int DarcyStokes::init(const MeshInfo& mi, PhysProperty * pp, const std::vector<double>& param){

    basis_ = basis();
    hdiv_  = Hdivmixed();
    br_    = BRMixed();

    br_.ComputeTotalDOF(mi);
    hdiv_.ComputeTotalDOF(mi);

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    totalElem = M * N;

    int tolgauss = ((M+1)*N + M*(N+1))*3;

    StokesVel.resize(tolgauss);
    DarcyVel.resize(tolgauss);

    // Compute boundary conditions
    ComputeEssenBndryAll(mi,pp,param);
    ComputeNaturBndryAll(mi,pp,param);

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

    cout << "Matrix Kg : " << endl;
    MatView(redsys.Kg, PETSC_VIEWER_STDOUT_WORLD);

    cout << "Matrix B : " << endl;
    MatView(redsys.B, PETSC_VIEWER_STDOUT_WORLD);

    cout << "Matrix Bg : " << endl;
    MatView(redsys.Bg, PETSC_VIEWER_STDOUT_WORLD);

    cout << "Matrix C : " << endl;
    MatView(redsys.C, PETSC_VIEWER_STDOUT_WORLD);

    cout << "Essen boundary : " << endl;    
    VecView(redsys.g, PETSC_VIEWER_STDOUT_WORLD);

    cout << "Natural boundary : " << endl;
    VecView(redsys.neum, PETSC_VIEWER_STDOUT_WORLD);

    return 1;
}

int DarcyStokes::showMatrix(){

    printRedSys(reducedDarcy_); 

    return 1;
}

int DarcyStokes::ReconstructEdgeVel(const vector<vertex>& edgegaussp,
                                    const MeshInfo& mi){

    // Extract nest vectors
    Vec stokesv, darcyv;
    PetscCall(VecNestGetSubVec(result.x, 0, &stokesv));
    PetscCall(VecNestGetSubVec(result.x, 1, &darcyv));

    ExtractVelocityEdge(StokesVel, edgegaussp, mi, refArrayStokesEssen_ , 
                      &stokesv, &reducedStokes_.g, br_, basis_);

    ExtractVelocityEdge(DarcyVel, edgegaussp, mi, refArrayDarcyEssen_, 
						    &darcyv, &reducedDarcy_.g, hdiv_, basis_);

    return 1;
}
