#include <iostream>
#include <petsc.h>
#include "integral.h"
#include "input.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "assemble.h"
#include "util.h"
#include "myFunc.h"
#include "bndry.h"
#include "solve.h"
#include "error.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

int main(int argc, char ** argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Testing coupled Stokes-Darcy system \n"));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));

    // =================================================================================

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 2, N = 2;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    // Create data management object
    DM    dm;
    Vec   fullmesh;
    const int stencilWidth = 1;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidth, NULL, NULL, &dm);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dm);               CHKERRQ(ierr);
    ierr = DMSetUp(dm);                        CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm, &fullmesh);CHKERRQ(ierr); 

    PhysProperty * physproperty = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(physproperty);

    double physscale = physproperty->L0/physproperty->l0;
    double L = 2*physscale, H = 1*physscale;
    double xstart = -1*physscale, ystart = -1.0001*physscale;
//    double L = 2, H = 1;
//    double xstart = -1, ystart = -1.1;

    ierr = PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL); CHKERRQ(ierr);

    int singleStencilTest = 0;
    double scale = 1;
    ierr = PetscOptionsGetInt(NULL,NULL, "-single", &singleStencilTest, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL, "-scale", &scale, NULL);CHKERRQ(ierr);

    if (singleStencilTest){
        L = L/scale;
        H = H/scale;
        xstart = -L/2.0;
        ystart = -H/2.0;
    }

    MeshParam mp;
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    // Uniform or distorted mesh
    int meshtype=0;
    ierr = PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshtype,NULL);CHKERRQ(ierr);
    switch(meshtype){
        case 0: CreateFullMesh(dm, &fullmesh, &mp); break;
        case 1: LogicRectMesh(dm, &fullmesh, &mp);  break;
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    int printmesh=0;
    ierr = PetscOptionsGetInt(NULL,NULL,"-printmesh",&printmesh,NULL);CHKERRQ(ierr);
    if(printmesh){ 
        VecView(fullmesh, PETSC_VIEWER_STDOUT_WORLD);
        PrintFullMesh(dm, &fullmesh);
    }

    //cout << "Mesh Created. To check full mesh, rerun with -printmesh 1 " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // =================================================================================

    // Contain defined mesh in vector container.
    // and verify it.
    vector< valarray<double> > mesh;
    
    ReadMeshPortion(dm, &fullmesh, mesh);

    //cout << "Converted c array of local mesh into vector container c++ " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // =================================================================================

    DM dmu;

    int cell_ghost = 0;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, cell_ghost, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu);               CHKERRQ(ierr);
    ierr = DMSetUp(dmu);                        CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    Vec localu; 
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu);

    // It can be changed later to not be double
    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    // =================================================================================
    // Create MeshInfo object
    MeshInfo mi; 

    // Assign local mesh and local values to mi
    mi.lmesh = mesh;
    mi.localVals = lu;

    AssignValuesMeshInfo(mi,dm,dmu); 

    // =================================================================================
    // Declare classes for shape functions
    basis * basis_   = new basis();
    Hdivmixed * hdiv = new Hdivmixed();
    BRMixed * br     = new BRMixed();

    br->ComputeTotalDOF(mi);
    hdiv->ComputeTotalDOF(mi);

    // Assemble matrices, create full system
    System * system = (System *)malloc(sizeof(System));

    SerialMatrixAssembleBlock(mi, *basis_, *hdiv, *br, physproperty, system);

    // Assign boundary conditions
    bndryVal bndryStokesDiri;
    bndryVal bndryStokesNeum;
    bndryVal bndryDarcyDiri;
    bndryVal bndryDarcyNeum;

    MarkBndryDOFStokes(bndryStokesDiri, bndryStokesNeum, mi, *basis_, *br, physproperty);
    MarkBndryDOFDarcy(bndryDarcyDiri, bndryDarcyNeum, mi, *basis_, *hdiv, physproperty); 

    // Create reduced system
    ReducedSys * reducedDarcy = (ReducedSys *)malloc(sizeof(ReducedSys));
    ReducedSys * reducedStokes = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateReducedSerial(reducedDarcy, &system->Ad, &system->Bd, &system->sourceDarcy, bndryDarcyDiri);
    CreateReducedSerial(reducedStokes, &system->As, &system->Bs, &system->sourceStokes, bndryStokesDiri);

    CreateNeumBndryVec(br->getDOF(), bndryStokesDiri.size(), reducedStokes, bndryStokesNeum, bndryStokesDiri);

    CreateNeumBndryVec(hdiv->getDOF(), bndryDarcyDiri.size(), reducedDarcy, bndryDarcyNeum, bndryDarcyDiri);

    PetscCall(VecAXPY(reducedStokes->source, 1.0, reducedStokes->neum));

    linearSys * lsStokes = (linearSys *)malloc(sizeof(linearSys));
    linearSys * lsDarcy  = (linearSys *)malloc(sizeof(linearSys));

    CreateLinearSys(lsStokes, reducedStokes);
    CreateLinearSys(lsDarcy, reducedDarcy);

    // Check linear system component
    const char *checkAd = "MatrixCheckAd.dat";
    // Write A matrix
    WriteMat(lsDarcy->A,checkAd);

    const char *checkBd = "MatrixCheckBd.dat";
    // Write B matrix
    WriteMat(lsDarcy->B,checkBd);

    // Write right two right hand side vectors
    const char *checkgd1 = "MatrixCheckgd1.dat";
    WriteVec(lsDarcy->f, checkgd1);   

    const char *checkgd2 = "MatrixCheckgd2.dat";
    WriteVec(lsDarcy->g, checkgd2);   

    const char *checkAs = "MatrixCheckAs.dat";
    // Write A matrix
    WriteMat(lsStokes->A,checkAs);

    const char *checkBs = "MatrixCheckBs.dat";
    // Write B matrix
    WriteMat(lsStokes->B,checkBs);

    // Write right two right hand side vectors
    const char *checkgs1 = "MatrixCheckgs1.dat";
    WriteVec(lsStokes->f, checkgs1);   

    const char *checkgs2 = "MatrixCheckgs2.dat";
    WriteVec(lsStokes->g, checkgs2);   

//    VecView(lsStokes->g, PETSC_VIEWER_STDOUT_WORLD);
    // Solve a coupled system ==========================================================
    // Couple two saddle point system
    // Control number of iterations and tolerance
    int maxIter = 1;
    PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL);

    double tauUzawa1 = 10;
    double tauUzawa2 = 1;
    PetscOptionsGetReal(NULL, NULL, "-tau1", &tauUzawa1, NULL);
    PetscOptionsGetReal(NULL, NULL, "-tau2", &tauUzawa2, NULL);

    double tolUzawa = 10e-7;
    PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL);

    // Assign correct C matrix to the target system
    PetscCall(MatConvert(system->Cd, MATSAME, MAT_INITIAL_MATRIX, &lsDarcy->C));
    PetscCall(MatConvert(system->Cs, MATSAME, MAT_INITIAL_MATRIX, &lsStokes->C));

    //MatView(system->Cd, PETSC_VIEWER_STDOUT_WORLD);

    const char *checkCd = "MatrixCheckCd.dat";
    WriteMat(system->Cd,checkCd);

    const char *checkCs = "MatrixCheckCs.dat";
    WriteMat(system->Cs,checkCs);

    const char *checkK = "MatCheckK.dat";
    WriteMat(system->K,checkK);

    PetscCall(VecZeroEntries(lsDarcy->x));
    PetscCall(VecZeroEntries(lsDarcy->y));
    PetscCall(VecZeroEntries(lsStokes->x));
    PetscCall(VecZeroEntries(lsStokes->y));

    linearSys * lsResult = (linearSys *)malloc(sizeof(linearSys));

    int precondType = 1;
    PetscOptionsGetInt(NULL, NULL, "-pType", &precondType, NULL);

    // ================================================================
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Mesh size : %d, %d.\nUzawa tolerance : %f \nMaximum Iteration : %d \nPreconditioner type : %d \n", M, N, tolUzawa,maxIter, precondType));
    switch (precondType){
        case 0:
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Single tau : %f\n", tauUzawa1));
        break;

        case 1:
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Split tau : %f, %f\n", tauUzawa1, tauUzawa2));
        break;

        case 2:
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Stokes tau : %f\n", tauUzawa1));
        break;

        default:
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Incorrect preconditioner type %d \n",precondType));
        break;
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));
    // ================================================================

    CoupledSolver(lsStokes, lsDarcy, lsResult, &system->K, tolUzawa, maxIter, tauUzawa1, tauUzawa2, precondType);
    //CoupledSolver(lsStokes, lsDarcy, lsResult, &system->K, tolUzawa, maxIter, tauUzawa1);


    // Result output - Stokes and Darcy
    Vec stokesv;
    Vec darcyv;

    PetscCall(VecNestGetSubVec(lsResult->x, 0, &stokesv));
    PetscCall(VecNestGetSubVec(lsResult->x, 1, &darcyv));

    std::vector<double> fullsolStokes = GetFullSol(&stokesv, bndryStokesDiri, br->getDOF());
    std::vector<double> fullsolDarcy  = GetFullSol(&darcyv, bndryDarcyDiri, hdiv->getDOF());

    // Stokes quiver output
    quiverOutput(mi, fullsolStokes, fullsolDarcy, M, N, *basis_, *br, *hdiv, physproperty);

    // =============================================================================
    // Result output - Single Stokes
    //std::vector<double> fullsolStokes = GetFullSol(&lsResult->x, bndryStokesDiri, br->getDOF());

    //quiverOutput(mi, fullsolStokes, M, N, *basis_,*br, physproperty);

    // =================================================================================
    // Clear used objects
    DMDAVecRestoreArray(dmu,localu,&lu);
    DMRestoreLocalVector(dmu, &localu); 

    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    DMDestroy(&dm);
    DMDestroy(&dmu);

    PetscFinalize();

    return 0;

}
