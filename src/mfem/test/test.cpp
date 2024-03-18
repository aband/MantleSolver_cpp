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

    // =================================================================================

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

    PhysProperty * physproperty = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(physproperty);

//    double L = 2*160000/physproperty->x0, H = 1*160000/physproperty->x0;
//    double xstart = -1*160000/physproperty->x0, ystart = 0;
    double L = 2, H = 1;
    double xstart = -1, ystart = -1.1;
//    double L = 2.0, H = 2.0;
//    double xstart = -1.0, ystart = -1.0;

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

    // Define basis functions H(div) conforming and Direct Serendipity space
    basis * testBasis = new basis();

    Hdivmixed * hdiv = new Hdivmixed();

    BRMixed * br = new BRMixed();

    System * system = (System *)malloc(sizeof(System));

    // Allocate space for matrix struct

    br->ComputeTotalDOF(mi);

    hdiv->ComputeTotalDOF(mi);

    // Create two physics system at the same time
    // Darcy and Stokes systems
    SerialMatrixAssembleBlock(mi, *testBasis, *hdiv, *br, physproperty, system);

    // Write coupling and compaction matrix ============================================
    const char *checkCd = "MatrixCheckCd.dat";
    WriteMat(system->Cd,checkCd);

    const char *checkCs = "MatrixCheckCs.dat";
    WriteMat(system->Cs,checkCs);

    const char *checkK = "MatCheckK.dat";
    WriteMat(system->K,checkK);
    // =================================================================================

    // Create reduced system
    // right hand side vectors stem from the created reduced system 
    // Reduced system for Stokes and Darcy sytem are being created separately
    bndryVal bndryStokes;
    bndryVal bndryDarcy;

    MarkBndryDOFDarcy(bndryDarcy, mi, (*testBasis), (*hdiv), physproperty);
    MarkBndryDOFStokes(bndryStokes, mi, (*testBasis), (*br), physproperty);

    // Separate Dirichlet and Neumann boundary condition
    bndryVal bndryStokesDiri;
    bndryVal bndryStokesNeum;

    MarkBndryDOFStokes(bndryStokesDiri,bndryStokesNeum,mi,(*testBasis),(*br),physproperty);

    // Test reduced system
    ReducedSys * redTest = (ReducedSys *)malloc(sizeof(ReducedSys));
    CreateReducedSerial(redTest, &system->As, &system->Bs, &system->sourceStokes, bndryStokesDiri);

    // ==================================================================================

    ReducedSys * reducedsys = (ReducedSys *)malloc(sizeof(ReducedSys));

    // Two systems are created at the same time.
    // returns A, B, C, K matrices
    CreateReducedSerial(reducedsys, &system->Ad, &system->Bd, &system->sourceDarcy, bndryDarcy);

    linearSys * ls = (linearSys *)malloc(sizeof(linearSys));

    // Create Target linear system
    CreateLinearSys(ls, reducedsys);

    //PetscCall(MatConvert(system->Cd, MATSAME, MAT_INITIAL_MATRIX, &ls->C));

    // Test inexect Uzawa iteration algorithm
    //InexactUzawa(ls, 10e-10, 30000, 0.08);

    // Write generated reduced matrices ================================================
    const char *checkAd = "MatrixCheckAd.dat";
    // Write A matrix
    WriteMat(ls->A,checkAd);

    const char *checkBd = "MatrixCheckBd.dat";
    // Write B matrix
    WriteMat(ls->B,checkBd);

    // Write right two right hand side vectors
    const char *checkgd1 = "MatrixCheckgd1.dat";
    WriteVec(ls->f, checkgd1);   

    const char *checkgd2 = "MatrixCheckgd2.dat";
    WriteVec(ls->g, checkgd2);   
    // =================================================================================
/*
    // Write out L2 error for Darcy only system
    std::vector<double> fullSolDarcyOnly = GetFullSol(&ls->x,bndryDarcy,hdiv->getDOF());
    double errorSumuDarcyOnly = 0.0;
        for (int j=0; j<N; j++){
        for (int i=0; i<M; i++){

            testBasis->GetCorners(mi,{i,j});

            // Darcy
            std::array<double, 8> sewD = ExtractWeights(fullSolDarcyOnly, hdiv->LocalToGlobal(mi,{i,j})); 
            errorSumuDarcyOnly += L2ErrorElem(sewD, {i,j}, bndryu,physproperty, gwf, gpf, *testBasis, *hdiv);
}}
        cout << "Darcy Only: ||u-u_h||_L2 : " <<  pow(errorSumuDarcyOnly,0.5) << endl;
*/
    // =================================================================================
    // Test of stokes equation starts from here
    ReducedSys * reducedsysStokes = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateReducedSerial(reducedsysStokes, &system->As, &system->Bs, &system->sourceStokes, bndryStokes);

    linearSys * lsStokes = (linearSys *)malloc(sizeof(linearSys));

    CreateLinearSys(lsStokes, reducedsysStokes);

    //PreconditionedUzawa(lsStokes, 10e-10, 30000, 1);

    // Write out linear system =========================================================
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
    // =================================================================================
/*
    std::vector<double> fullSolStokesOnly = GetFullSol(&lsStokes->x, bndryStokes, br->getDOF());

    double errorSumuStokesOnly = 0.0;

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        testBasis->GetCorners(mi,{i,j});

        // Stokes
        std::array<double, 12> sewStokes = ExtractWeights(fullSolStokesOnly, br->LocalToGlobal(mi,{i,j})); 
        errorSumuStokesOnly += L2ErrorElem(sewStokes, {i,j}, bndryVs, physproperty, gwf, gpf, *testBasis, *br);

}}
    cout << "Stokes  Only: ||u-u_h||_L2 : " <<  pow(errorSumuStokesOnly,0.5) << endl;
*/

    // Solve a coupled system ==========================================================
    // Couple two saddle point system
    // Control number of iterations and tolerance
    int maxIter;
    PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL);

    double tauUzawa;
    PetscOptionsGetReal(NULL, NULL, "-tau", &tauUzawa, NULL);

    double tolUzawa;
    PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL);

    // Assign correct C matrix to the target system
    PetscCall(MatConvert(system->Cd, MATSAME, MAT_INITIAL_MATRIX, &ls->C));
    PetscCall(MatConvert(system->Cs, MATSAME, MAT_INITIAL_MATRIX, &lsStokes->C));

    //VecView(ls->x,PETSC_VIEWER_STDOUT_WORLD);
    //VecView(ls->y,PETSC_VIEWER_STDOUT_WORLD);

    VecZeroEntries(ls->x);
    VecZeroEntries(ls->y);
    VecZeroEntries(lsStokes->x);
    VecZeroEntries(lsStokes->y);

    linearSys * lsResult = (linearSys *)malloc(sizeof(linearSys));

    CoupledSolver(lsStokes, ls, lsResult, &system->K, tolUzawa, maxIter, tauUzawa);

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
    int quiverType = 1;
    PetscOptionsGetInt(NULL, NULL, "-qt", &quiverType, NULL);

    if (checkError){

        // extract sub vectors from nest vector
        Vec stokesx;
        Vec darcyx;

        VecNestGetSubVec(lsResult->x, 0, &stokesx);
        VecNestGetSubVec(lsResult->x, 1, &darcyx);

        std::vector<double> fullSolStokes;
        std::vector<double> fullSolDarcy;
        //fullSol = GetFullSol(&ls->x,bndryDarcy,hdiv->getDOF());
        //fullSol = GetFullSol(&testReduced, bndryTest, hdiv->getDOF());
        //fullSol = GetFullSol(&lsStokes->x, bndryStokes, br->getDOF());

        fullSolStokes = GetFullSol(&stokesx, bndryStokes, br->getDOF());
        fullSolDarcy  = GetFullSol(&darcyx, bndryDarcy, hdiv->getDOF());

        // Fetch gauss points and gauss weights
        double errorSumu = 0.0;
        double errorSump = 0.0;

        double errorSumuStokes = 0.0;
        double errorSumuDarcy  = 0.0;

        //double *arrayp;
        //PetscCall(VecGetArray(lsStokes->y,&arrayp));
        //PetscCall(VecGetArray(ls->y,&arrayp));

        for (int j=0; j<N; j++){
        for (int i=0; i<M; i++){

            testBasis->GetCorners(mi,{i,j});

            // Coupled
            std::array<double, 12> singleWgtsStokes = ExtractWeights(fullSolStokes, br->LocalToGlobal(mi,{i,j}));
            std::array<double, 8> singleWgtsDarcy = ExtractWeights(fullSolDarcy, hdiv->LocalToGlobal(mi,{i,j}));

            errorSumuStokes += L2ErrorElem(singleWgtsStokes,{i,j},bndryVs,physproperty,gwf,gpf,*testBasis,*br);
            errorSumuDarcy  += L2ErrorElem(singleWgtsDarcy, {i,j},bndryu, physproperty,gwf,gpf,*testBasis,*hdiv);

            // Output approximation and exact velocity for quiver plot

            // Center of element
            //cout << "( " << j << ", " << i << ") : " << errorSumuStokes  << ", " << errorSumuDarcy << " ";
}}
//        }cout << endl; }
        //PetscCall(VecRestoreArray(lsStokes->y,&arrayp));
        //PetscCall(VecRestoreArray(ls->y,&arrayp));

        //cout << "||u-u_h||_L2 : " <<  pow(errorSumu,0.5) << endl;
        //cout << "||p-p_h||_L2 : " <<  pow(errorSump,0.5) << endl;
        cout << "Coupled : ||u-h_h||_L2 : " << pow(errorSumuStokes+errorSumuDarcy,0.5) << endl;
        cout << "Stokes  : ||u-u_h||_L2 : " <<  pow(errorSumuStokes,0.5) << endl;
        cout << "Darcy   : ||u-u_h||_L2 : " <<  pow(errorSumuDarcy,0.5) << endl;

        if (quiverType == 1){
            quiverOutput(mi,fullSolDarcy,M,N,*testBasis,*br,*hdiv,physproperty,1);
        } else {
            quiverOutput(mi,fullSolStokes,M,N,*testBasis,*br,*hdiv,physproperty,2);
        }
    }

    Vec stokesp;
    Vec darcyp;

    VecNestGetSubVec(lsResult->y, 0, &stokesp);
    VecNestGetSubVec(lsResult->y, 1, &darcyp);

    double stokesmean = 0.0;
    double darcymean = 0.0;

    VecMean(stokesp, &stokesmean);
    VecMean(darcyp, &darcymean);

    cout << stokesmean << " " << darcymean << endl;
    //VecView(lsResult->y,PETSC_VIEWER_STDOUT_WORLD);
    //VecView(lsResult->x,PETSC_VIEWER_STDOUT_WORLD);

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
