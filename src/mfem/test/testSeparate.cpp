// Test Darcy ans Stokes system separately
// A decoupled system
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

int main(int argc, char **argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 2, N = 2;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size);CHKERRQ(ierr);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Mesh size : %d, %d.\n", M, N));

    ierr = PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n");CHKERRQ(ierr);

    // =================================================================================

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

    //double physscale = physproperty->L0/physproperty->l0;
    double physscale = 1.0;
    double L = 2*physscale, H = 1*physscale;
    double xstart = -1*physscale, ystart = -1.00*physscale;
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

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;

    // Define basis functions H(div) conforming space
    basis * testBasis = new basis();
    Hdivmixed * hdiv  = new Hdivmixed();
    BRMixed * br      = new BRMixed();

    // Allocate space for matrix struct
    br->ComputeTotalDOF(mi);

    hdiv->ComputeTotalDOF(mi);

    // Allocate memory space for linear system
    System * system = (System *)malloc(sizeof(System));

    SerialMatrixAssembleBlock(mi, *testBasis, *hdiv, *br, physproperty, system);

    // Mark boundary condition
    bndryVal bndryStokesDiri;
    bndryVal bndryStokesNeum;
    bndryVal bndryDarcyDiri;
    bndryVal bndryDarcyNeum;

    MarkBndryDOFStokes(bndryStokesDiri, bndryStokesNeum, mi, *testBasis, *br, physproperty);
    MarkBndryDOFDarcy(bndryDarcyDiri, bndryDarcyNeum, mi, *testBasis, *hdiv, physproperty); 

    // Create reduced system

    ReducedSys * reducedDarcy = (ReducedSys *)malloc(sizeof(ReducedSys));
    ReducedSys * reducedStokes = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateReducedSerial(reducedDarcy, &system->Ad, &system->Bd, &system->sourceDarcy, bndryDarcyDiri);
    CreateReducedSerial(reducedStokes, &system->As, &system->Bs, &system->sourceStokes, bndryStokesDiri);

    CreateNeumBndryVec(br->getDOF(), bndryStokesDiri.size(), reducedStokes, bndryStokesNeum, bndryStokesDiri);

    CreateNeumBndryVec(hdiv->getDOF(), bndryDarcyDiri.size(), reducedDarcy, bndryDarcyNeum, bndryDarcyDiri);

    // Test Dirichlet problem first
    //PetscCall(VecAXPY(reducedStokes->source, 1.0, reducedStokes->neum));

    linearSys * lsStokes = (linearSys *)malloc(sizeof(linearSys));
    linearSys * lsDarcy  = (linearSys *)malloc(sizeof(linearSys));

    CreateLinearSys(lsStokes, reducedStokes);
    CreateLinearSys(lsDarcy, reducedDarcy);

    PetscCall(MatConvert(system->Cd, MATSAME, MAT_INITIAL_MATRIX, &lsDarcy->C));
    PetscCall(MatConvert(system->Cs, MATSAME, MAT_INITIAL_MATRIX, &lsStokes->C));

    PetscCall(MatZeroEntries(lsDarcy->C));
    PetscCall(MatZeroEntries(lsStokes->C));

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

    // =======================================================================
    int maxIterStokes = 1;
    PetscOptionsGetInt(NULL, NULL, "-maxIterStokes", &maxIterStokes, NULL);

    int maxIterDarcy = 1;
    PetscOptionsGetInt(NULL, NULL, "-maxIterDarcy", &maxIterDarcy, NULL);

    double tauUzawaStokes = 20;
    PetscOptionsGetReal(NULL, NULL, "-tauStokes", &tauUzawaStokes, NULL);

    double tauUzawaDarcy = 1;
    PetscOptionsGetReal(NULL, NULL, "-tauDarcy", &tauUzawaDarcy, NULL);

    double tolUzawaStokes = 10e-7;
    PetscOptionsGetReal(NULL, NULL, "-tolStokes", &tolUzawaStokes, NULL);

    double tolUzawaDarcy = 10e-7;
    PetscOptionsGetReal(NULL, NULL, "-tolDarcy", &tolUzawaDarcy, NULL);

    int precondType = 0;
    PetscOptionsGetInt(NULL, NULL, "-pType", &precondType, NULL);

    // Solve Stokes and Darcy sub problems separately
    // Using simple uzawa here with C = 0 (no coupling matrix)

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Summary of Darcy problem. \nUzawa tolerance : %f \nMaximum Iteration : %d \nTau : %f\n", tolUzawaDarcy,maxIterDarcy,tauUzawaDarcy));
    //SimpleUzawa(lsDarcy, tolUzawaDarcy, maxIterDarcy, tauUzawaDarcy, precondType);
    ExactUzawa(lsDarcy, tolUzawaDarcy, maxIterDarcy);

    // ============================================================================
    // Write out L2 error for Darcy only system
    std::vector<double> fullSolDarcy = GetFullSol(&lsDarcy->x,bndryDarcyDiri,hdiv->getDOF());
    double errorSumuDarcy = 0.0;
        for (int j=0; j<N; j++){
        for (int i=0; i<M; i++){

            testBasis->GetCorners(mi,{i,j});

            // Darcy
            std::array<double, 8> sewD = ExtractWeights(fullSolDarcy, hdiv->LocalToGlobal(mi,{i,j})); 
            errorSumuDarcy += L2ErrorElem(sewD, {i,j}, bndryu,physproperty, gwf, gpf, *testBasis, *hdiv);
}}

    cout << "Darcy Only: ||u-u_h||_L2 : " <<  pow(errorSumuDarcy,0.5) << endl;

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));

    // Write out L2 error for Stokes only system

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Summary of Stokes problem. \nUzawa tolerance : %f \nMaximum Iteration : %d \nTau : %f\n", tolUzawaStokes,maxIterStokes,tauUzawaStokes));

    SimpleUzawa(lsStokes, tolUzawaStokes, maxIterStokes, tauUzawaStokes, 0);

    std::vector<double> fullSolStokes = GetFullSol(&lsStokes->x, bndryStokesDiri, br->getDOF());

    double errorSumuStokes = 0.0;

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        testBasis->GetCorners(mi,{i,j});

        // Stokes
        std::array<double, 12> sewStokes = ExtractWeights(fullSolStokes, br->LocalToGlobal(mi,{i,j})); 
        errorSumuStokes += L2ErrorElem(sewStokes, {i,j}, bndryVs, physproperty, gwf, gpf, *testBasis, *br);

}}

    cout << "Stokes Only: ||u-u_h||_L2 : " <<  pow(errorSumuStokes,0.5) << endl;

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));

    return 0;
}
