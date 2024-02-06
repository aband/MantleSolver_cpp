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

    ierr = PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n");CHKERRQ(ierr);

    // ==========================================================================================================================

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 3, N = 3;
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

    double L = 2.0, H = 2.0;
    double xstart = -1.0, ystart = -1.0;
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

    // ==========================================================================================================================

    // Contain defined mesh in vector container.
    // and verify it.
    vector< valarray<double> > mesh;
    
    ReadMeshPortion(dm, &fullmesh, mesh);

    //cout << "Converted c array of local mesh into vector container c++ " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // ====================================================================================================================================

    DM dmu;

    int cell_ghost = 0;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, cell_ghost, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu);               CHKERRQ(ierr);
    ierr = DMSetUp(dmu);                        CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    // Initialize with oblique data for Burgers equation 
    //ObliqueBurgers(dm,dmu,&fullmesh,&globalu,Initial_Condition);
    //SimpleInitialValue(dm,dmu,&fullmesh,&globalu,{-L/(2*M)},func);

    Vec localu; 
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu);

    // It can be changed later to not be double
    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    // ====================================================================================================================================

    // Create MeshInfo object
    MeshInfo mi; 

    // Assign local mesh and local values to mi
    mi.lmesh = mesh;
    mi.localVals = lu;

    AssignValuesMeshInfo(mi,dm,dmu); 

// ========================================================================================================================================

    // Define basis functions H(div) conforming and Direct Serendipity space
    basis * testBasis = new basis();

    Hdivmixed * hdiv = new Hdivmixed();

    BRMixed * br = new BRMixed();

    System * system = (System *)malloc(sizeof(System));

    PhysProperty * physproperty = (PhysProperty *)malloc(sizeof(PhysProperty));

    (*physproperty).l = 20;

    // Allocate space for matrix struct

    br->ComputeTotalDOF(mi);

    hdiv->ComputeTotalDOF(mi);

    // Create two physics system at the same time
    // Darcy and Stokes systems
    SerialMatrixAssembleBlock(mi, *testBasis, *hdiv, *br, physproperty, system);

    // Create reduced system
    // right hand side vectors stem from the created reduced system 
    // Reduced system for Stokes and Darcy sytem are being created separately
    bndryVal bndryStokes;
    bndryVal bndryDarcy;

    MarkBndryDOFDarcy(bndryDarcy, mi, (*testBasis), (*hdiv));
    MarkBndryDOFStokes(bndryStokes, mi, (*testBasis), (*br));

    ReducedSys * reducedsys = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateReducedSerial(reducedsys, &system->Ad, &system->Bd, &system->sourceDarcy, bndryDarcy);

    // Test Darcy part alone ============================================================
    Vec g1;

    // Check mat size
    int cM, cN;
    PetscCall(MatGetSize(reducedsys->M, &cM, &cN));

    PetscCall(VecCreate(PETSC_COMM_WORLD, &g1));
    PetscCall(VecSetSizes(g1, PETSC_DECIDE, cM));
    PetscCall(VecSetUp(g1));
  
    PetscCall(MatMult(reducedsys->Kg, reducedsys->g, g1));

    Vec g2;
    PetscCall(MatGetSize(reducedsys->B, &cM, &cN));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &g2));
    PetscCall(VecSetSizes(g2, PETSC_DECIDE, cN));
    PetscCall(VecSetUp(g2));

    Mat BgT;
    PetscCall(MatCreateTranspose(reducedsys->Bg, &BgT));
 
    PetscCall(MatMult(BgT, reducedsys->g, g2));

    // Move boundary condition vectors to the right hand side of the 
    //VecScale(g1,-1);
    PetscCall(VecAYPX(g1,-1,reducedsys->source));
    PetscCall(VecScale(g2,-1));

// Check computed system
/*
    const char *checkA = "MatrixCheckA.dat";
    // Write A matrix
    WriteMat(reducedsys->M,checkA);

    const char *checkB = "MatrixCheckB.dat";
    // Write B matrix
    WriteMat(reducedsys->B,checkB);

    // Write right two right hand side vectors
    const char *checkg1 = "MatrixCheckg1.dat";
    WriteVec(g1, checkg1);   

    const char *checkg2 = "MatrixCheckg2.dat";
    WriteVec(g2, checkg2);   

    const char *checkKg = "MatrixCheckKg.dat";
    WriteMat(reducedsys->Kg,checkKg);

    const char *checkg = "VecCheckg.dat";
    WriteVec(reducedsys->g,checkg);
*/
    // Test inexect Uzawa iteration algorithm
    linearSys * ls = (linearSys *)malloc(sizeof(linearSys));

    PetscCall(MatConvert(reducedsys->B, MATSAME, MAT_INITIAL_MATRIX, &ls->B));
    PetscCall(MatConvert(reducedsys->M, MATSAME, MAT_INITIAL_MATRIX, &ls->A));

    PetscCall(VecDuplicate(g1,&ls->f));
    PetscCall(VecDuplicate(g2,&ls->g));

    VecCopy(g1, ls->f);
    VecCopy(g2, ls->g);

    PetscCall(VecCreate(PETSC_COMM_WORLD, &ls->x));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &ls->y));

    PetscCall(VecSetSizes(ls->x,PETSC_DECIDE,cM));
    PetscCall(VecSetSizes(ls->y,PETSC_DECIDE,cN));

    PetscCall(VecSetUp(ls->x));
    PetscCall(VecSetUp(ls->y));

    PetscCall(VecCopy(g1,ls->x));
    PetscCall(VecCopy(g2,ls->y));
   
    PetscCall(VecZeroEntries(ls->x));
    PetscCall(VecZeroEntries(ls->y));

    // Control number of iterations and tolerance
    int maxIter;
    PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL);

    double tauUzawa;
    PetscOptionsGetReal(NULL, NULL, "-tau", &tauUzawa, NULL);

    if (tauUzawa < 0){
        // Use element size related tauUzawa
        tauUzawa = 1.0/(double)N / (double) M;
    }

    double tolUzawa;
    PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL);

    //PreconditionedUzawa(ls, tolUzawa, maxIter, tauUzawa);

    // =================================================================================
    // End of test of Darcy equation (literally Poisson equation 
    // turns second order equation into first linear system )
    // Test of stokes equation starts from here
    ReducedSys * reducedsysStokes = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateReducedSerial(reducedsysStokes, &system->As, &system->Bs, &system->sourceStokes, bndryStokes);

//    VecView(system->sourceStokes, PETSC_VIEWER_STDOUT_WORLD);

    Vec g1Stokes;
    PetscCall(MatGetSize(reducedsysStokes->M, &cM, &cN));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &g1Stokes));
    PetscCall(VecSetSizes(g1Stokes,PETSC_DECIDE,cM));
    PetscCall(VecSetUp(g1Stokes));

    PetscCall(MatMult(reducedsysStokes->Kg, reducedsysStokes->g, g1Stokes));
 
    Vec g2Stokes;
    PetscCall(MatGetSize(reducedsysStokes->B, &cM, &cN));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &g2Stokes));
    PetscCall(VecSetSizes(g2Stokes,PETSC_DECIDE,cN));
    PetscCall(VecSetUp(g2Stokes));

    Mat BgTStokes;
    PetscCall(MatCreateTranspose(reducedsysStokes->Bg, &BgTStokes));
    PetscCall(MatMult(BgTStokes, reducedsysStokes->g, g2Stokes));

    PetscCall(VecAYPX(g1Stokes, -1, reducedsysStokes->source));
    PetscCall(VecScale(g2Stokes, -1));

    linearSys * lsStokes = (linearSys *)malloc(sizeof(linearSys));

    PetscCall(MatConvert(reducedsysStokes->B, MATSAME, MAT_INITIAL_MATRIX, &lsStokes->B));
    PetscCall(MatConvert(reducedsysStokes->M, MATSAME, MAT_INITIAL_MATRIX, &lsStokes->A));

    PetscCall(VecDuplicate(g1Stokes,&lsStokes->f));
    PetscCall(VecDuplicate(g2Stokes,&lsStokes->g));

    PetscCall(VecCopy(g1Stokes, lsStokes->f));
    PetscCall(VecCopy(g2Stokes, lsStokes->g));

    PetscCall(VecCreate(PETSC_COMM_WORLD, &lsStokes->x));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &lsStokes->y));

    PetscCall(VecSetSizes(lsStokes->x,PETSC_DECIDE,cM));
    PetscCall(VecSetSizes(lsStokes->y,PETSC_DECIDE,cN));

    PetscCall(VecSetUp(lsStokes->x));
    PetscCall(VecSetUp(lsStokes->y));

    PetscCall(VecCopy(g1Stokes,lsStokes->x));
    PetscCall(VecCopy(g2Stokes,lsStokes->y));
   
    PetscCall(VecZeroEntries(lsStokes->x));
    PetscCall(VecZeroEntries(lsStokes->y));

    //PetscCall(MatConvert(system->Cs, MATSAME, MAT_INITIAL_MATRIX, &lsStokes->C));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &lsStokes->C));
    PetscCall(MatSetSizes(lsStokes->C, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N));
    PetscCall(MatSetUp(lsStokes->C));

    PetscCall(MatZeroEntries(lsStokes->C));

    PreconditionedUzawa(lsStokes, tolUzawa, maxIter, tauUzawa);

    const char *checkA = "MatrixCheckA.dat";
    // Write A matrix
    WriteMat(reducedsysStokes->M,checkA);

    const char *checkB = "MatrixCheckB.dat";
    // Write B matrix
    WriteMat(reducedsysStokes->B,checkB);

    // Write right two right hand side vectors
    const char *checkg1 = "MatrixCheckg1.dat";
    WriteVec(g1Stokes, checkg1);   

    const char *checkg2 = "MatrixCheckg2.dat";
    WriteVec(g2Stokes, checkg2);   

    const char *checkKg = "MatrixCheckKg.dat";
    WriteMat(reducedsysStokes->Kg,checkKg);

    const char *checkg = "VecCheckg.dat";
    WriteVec(reducedsysStokes->g,checkg);

    // =================================================================================
    // Check solution created
/*    Vec testFull;
    Vec testReduced;

    double * arraytestfull;
    double * arraytestreduced;

    VecCreate(PETSC_COMM_WORLD, &testFull);
    VecCreate(PETSC_COMM_WORLD, &testReduced);

    VecSetSizes(testFull, PETSC_DECIDE, hdiv->getDOF());
    VecSetSizes(testReduced, PETSC_DECIDE, cM);
    VecSetUp(testFull);
    VecSetUp(testReduced);

    VecGetArray(testFull, &arraytestfull);
    VecGetArray(testReduced, &arraytestreduced);

    for (int i=0; i<hdiv->getDOF(); i++){
        arraytestfull[i] = i;
    }

    // Create manufactured boundary values and corresponding data structure
    bndryVal bndryTest;

    for (const auto& bv : bndryDarcy){
        bndryTest.insert(std::make_pair<int, bndryInfo>
                         ((int)bv.first, {0,(double)bv.first,{0,0}}));

    }

    // Create Test reduced vector
    int count = 0; 
    for (int i=0; i<hdiv->getDOF(); i++){
   
        auto ifFind = bndryTest.find(i);
        if (ifFind == bndryTest.end()){
            arraytestreduced[count] = i;
            count ++; 
        }
    }

    VecRestoreArray(testFull, &arraytestfull);
    VecRestoreArray(testReduced, &arraytestreduced);
*/

    // Check computed error results
    int checkError = 0;
    PetscOptionsGetInt(NULL, NULL, "-checkError", &checkError, NULL);

    if (checkError){

        std::vector<double> fullSol;
        //fullSol = GetFullSol(&ls->x,bndryDarcy,hdiv->getDOF());
        //fullSol = GetFullSol(&testReduced, bndryTest, hdiv->getDOF());
        fullSol = GetFullSol(&lsStokes->x, bndryStokes, br->getDOF());

        // Fetch gauss points and gauss weights
        const valarray<double>& gwf = GaussWeightsFace;
        const vector<vertex>& gpf = GaussPointsFace;

        double errorSumu = 0.0;
        double errorSump = 0.0;

        double *arrayp;
        PetscCall(VecGetArray(lsStokes->y,&arrayp));

        for (int j=0; j<N; j++){
        for (int i=0; i<M; i++){

            testBasis->GetCorners(mi,{i,j});

            // Stokes
            std::array<double, 12> singleElemWeights = ExtractWeights(fullSol, br->LocalToGlobal(mi,{i,j})); 
            errorSumu += L2ErrorElem(singleElemWeights, {i,j}, trueSol, gwf, gpf, *testBasis, *br);
            errorSump += L2ErrorElem(arrayp[j*M+i],trueSol,gwf,gpf,*testBasis,mi.cellArea.at(j*M+i));

            // Darcy
            //std::array<double, 8> singleElemWeights = ExtractWeights(fullSol, hdiv->LocalToGlobal(mi,{i,j})); 
            //errorSumu += L2ErrorElem(singleElemWeights, {i,j}, trueSol, gwf, gpf, *testBasis, *hdiv);
            //errorSump += L2ErrorElem(arrayp[j*M+i],trueSol,gwf,gpf,*testBasis,mi.cellArea.at(j*M+i));

        }}

        PetscCall(VecRestoreArray(lsStokes->y,&arrayp));

        cout << "||u-u_h||_L2 : " <<  errorSumu << endl;
        cout << "||p-p_h||_L2 : " <<  errorSump << endl;
    }

    // =================================================================================
    // Check FE function space
    // Check element {0,0}
    // (designed for single element case
    int seed = 5; 
    // assume L = H here
    double DX = L/(double)	M;

    double h = DX / (double) seed;

    int k = 0;
    int shift = 0;
    PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);
    PetscOptionsGetInt(NULL,NULL,"-s",&shift,NULL);

    testBasis->GetCorners(mi,{0,0});

    //std::array<double, 8> fakeweight = ExtractWeights(fullSol, hdiv->LocalToGlobal(mi, {0,0}));;

    std::array<double, 8> fakeweight = {0,0,0,0,1,1,1,1};

	 /*
    for (int j=seed; j>-1; j--){
    for (int i=0; i<seed + 1; i++){
        //std::array<vertex, 8> tmp = hdiv->ComputeHdivmixed(*testBasis, {xstart+i*h, ystart+j*h});
        std::array<std::array<double, 4>, 12> tmp = br->ComputeGradBRmixed(*testBasis, {xstart+i*h, ystart+j*h});
        //vertex tmp = testBasis->dR(k, {xstart+i*h, ystart+j*h});

        //std::array<vertex, 12> tmp = br->ComputeBRmixed(*testBasis, {xstart+i*h, ystart+j*h});

        cout << "(" << tmp[k+4*shift][0] << ", " << tmp[k+4*shift][1] << ", " << 
                       tmp[k+4*shift][2] << ", " << tmp[k+4*shift][3] <<  ")  ";

        //cout << "( " << tmp[0] << ", " << tmp[1] << " )" ;

        //vertex sum {0.0,0.0};
        //for (int g=0;g<8;g++){
        //    sum += fakeweight[g]*tmp[g]; 
        //}
        //cout << "(" << sum[0] << ", " << sum[1] << ")  ";
 
    }cout << endl;}

	 */
// ====================================================================================================================================
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
