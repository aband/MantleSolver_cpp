#include "stencilpolynomial.h"
#include "reconstruction.h"
#include "petsc.h"
#include "input.h"
#include "error.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

#include "temp.h"

int main(int argc, char ** argv){

    // Integrated test with petsc and mesh functions

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    int M = 10, N = 10;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    double L = 3, H = 1;
    double xstart = 0.0, ystart = 0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    int stencilWidthMesh = 5; // Ghost layer thickness for vertex
    int stencilWidthU = 3;    // Ghost layer thickness for cell

    DM dmu;
    DM dmMesh;

    int meshType = 0; 
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshType,NULL));

    double dt = 0.05*1.0/(double)M;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-dt", &dt, NULL));

    int Nt = 10;
    ierr = PetscOptionsGetInt(NULL,NULL,"-Nt",&Nt,NULL);CHKERRQ(ierr);

    Nt *= M;

    double CFL = dt/(1.0/(double)M);

    cout << "dt, dh = " << dt << " , " << 1.0/(double)M << ". " << "CFL number is : " << CFL << endl;

    M *= 3;

    // Create dmMesh
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, &dmMesh));
    PetscCall(DMSetFromOptions(dmMesh));              
    PetscCall(DMSetUp(dmMesh));

    // Create dmU
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 
    stencilWidthU, NULL, NULL, &dmu));
    PetscCall(DMSetFromOptions(dmu));              
    PetscCall(DMSetUp(dmu));     

    // Create MeshParam object (historical object one time use only)
	 MeshParam mp; 
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    Vec globalmesh;
    // Create global vector containing mesh
    PetscCall(DMCreateGlobalVector(dmMesh, &globalmesh));

    switch(meshType){
        case 0: CreateFullMesh(dmMesh, &globalmesh, &mp); break;
        case 1: LogicRectMesh(dmMesh, &globalmesh, &mp);  break;
        case 2: RefineMesh(dmMesh, &globalmesh, &mp);
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    MeshInfo mi;

    ReadMeshPortion(dmMesh, &globalmesh, mi.lmesh);

    AssignValuesMeshInfo(mi, dmMesh, dmu);

    cout << "check point 1" << endl;

    mi.L = L;
    mi.H = H;

    multilevel ml = multilevel();

cout << "Start here " << endl;
    ml.addLevel("(5,5)", {5,5}, mi);
cout << "(5,5) prepared." << endl;
    ml.addLevel("(3,3)", {3,3}, mi);
cout << "(3,3) prepared." << endl;
//    ml.addLevel("(2,2)", {2,2}, mi);
//cout << "(2,2) prepared." << endl;
cout << "End here " << endl;

    double h0 = sqrt((L*H)/(double)(M*N));

    // Define a usage for this multilevel weno
    mluse use = mluse();

    // Test for nonlinear weighting
    unordered_map<std::string, vector<indice>> interior;
//    interior.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-1,-1} }));
//    interior.insert(std::make_pair<std::string, vector<indice>>("(2,2)", { {-1,-1}, {0,-1}, {0,0}, {-1,0} }));
//    use.setmethod("interior", interior);
//    use.setbias("interior");

    interior.insert(std::make_pair<std::string, vector<indice>>("(5,5)", { {-2,-2} }));
    interior.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-2,-2}, {0,-2}, {-2,0}, {0,0} }));

    use.setmethod("interior", interior);
    use.setbias("interior");

    // add biased (3,3) stencil to side cells
    unordered_map<std::string, vector<indice>> side;
    //side.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-2,-1}, {0,-1} }));
    //side.insert(std::make_pair<std::string, vector<indice>>("(2,2)", { {-1,-1}, {0,-1}, {0,0}, {-1,0} }));
    side.insert(std::make_pair<std::string, vector<indice>>("(5,5)", {{-2,0}, {2,0} }));
    side.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-2,-2}, {0,-2}, {0,0}, {-2,0} }));

    use.setmethod("side", side);
    use.setbias("side");

    // =================================================================
    Vec globalvec;

    PetscCall(DMCreateGlobalVector(dmu, &globalvec));

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalvec, {0.0,0.0}, func);

    // =================================================================

    //simpleSSP2RK(dt, Nt, &globalvec, mi, ml, use, dmu, dmMesh);
	 cout << "Time stepping starts. " << endl;
//    simpleRK(dt, Nt, &globalvec, mi, ml, use, dmu, dmMesh);
//    simpleSSP2RK(dt, Nt, &globalvec, mi, ml, use, dmu, dmMesh);
    simpleSSP3RK(dt, Nt, &globalvec, mi, ml, use, dmu, dmMesh);
    reconPlot(mi, ml, use, 1, &globalvec, true, dmu, h0);
//    exactSol(mi,2.0, func, 1, true);

    // =================================================================

    PetscCall(VecDestroy(&globalvec));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    return 1;
}
