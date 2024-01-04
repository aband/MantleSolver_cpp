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

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

double func(const vertex& point, const vector<double>& param){
	 if (point[0]<param[0]){
//		  return point[0]*point[0]+point[1]*point[1];
	     return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1);
//        return sin(point[0] + point[1] + 0.1);
	 } else {
//		  return point[0]*point[0]*point[1]*point[1] + 1.0;
	     return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1) + 10;
//        return sin(point[0] + point[1] + 0.1) + 10;
	 }

    //return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1);
    //return sin(point[0] + point[1] + 0.1);
    //return point[0]*point[0] + point[1]*point[1];
    //return 0.5;
    //return point[0] + point[1];

}

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
    SimpleInitialValue(dm,dmu,&fullmesh,&globalu,{-L/(2*M)},func);

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

    Matrix * matrix = (Matrix *)malloc(sizeof(Matrix));

    PhysProperty * physproperty = (PhysProperty *)malloc(sizeof(PhysProperty));

    (*physproperty).l = 20;

    // Allocate space for matrix struct

    br->ComputeTotalDOF(mi);

    hdiv->ComputeTotalDOF(mi);

    Vec source;
    PetscCall(VecCreate(PETSC_COMM_WORLD, &source));
    PetscCall(VecSetSizes(source, PETSC_DECIDE, br->getDOF()));
    PetscCall(VecSetUp(source));

    SerialMatrixAssembleBlock(mi, *testBasis, *hdiv, *br, physproperty, matrix, &source);

    //const char *check1 = "MatrixCheck.dat";

    // Check matrix shape
    //DrawMat(matrix->As,check1);

    CreateSchurComplement(matrix, M*N, br->getDOF(), hdiv->getDOF());

    //const char *check2 = "schur.dat";

    //MatView(matrix->G, PETSC_VIEWER_STDOUT_WORLD);

    bndryVal bndryStokes;
    bndryVal bndryDarcy;
    // Create right hand side vector
    //MarkBndryDOFStokes(bndryStokes, mi, (*br));
    MarkBndryDOFDarcy(bndryDarcy, mi, (*testBasis), (*hdiv));
    MarkBndryDOFStokes(bndryStokes, mi, (*testBasis), (*br));

    ReducedSys * reducedsys = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateReducedSerial(reducedsys, &matrix->Ad, &matrix->Bd, bndryDarcy);

    // Test Darcy part alone
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

    // Create Schur complement
    Mat sub[4];
    Mat S1, S2, Sp1, Sp2;

    sub[0] = reducedsys->M;
    sub[1] = reducedsys->B;
    Mat BT;

    PetscCall(MatCreate(PETSC_COMM_WORLD, &BT));
    PetscCall(MatSetSizes(BT, PETSC_DECIDE, PETSC_DECIDE, cN, cM));
    PetscCall(MatSetUp(BT));

    sub[2] = BT;

    // Create zero matrix
    Mat Z;
    PetscCall(MatCreate(PETSC_COMM_WORLD, &Z));
    PetscCall(MatSetSizes(Z, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N));
    PetscCall(MatSetUp(Z));

    PetscCall(MatZeroEntries(Z));
    sub[3] = Z;

    Mat G;
    MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, sub, &G);

    VecScale(g1,-1);
    VecScale(g2,-1);

    Vec g[2];
    g[0] = g1;
    g[1] = g2;

//    Vec rhs;
//    VecCreateNest(PETSC_COMM_WORLD, 2, NULL, g, &rhs);

    PetscCall(MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY));

//    MatView(G, PETSC_VIEWER_STDOUT_WORLD);
//    VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);

    // Set up linear solver
    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, G, G);
    KSPSetType(ksp, KSPCG);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);

    Vec rhs;
    VecCreate(PETSC_COMM_WORLD, &rhs);
    VecSetSizes(rhs, PETSC_DECIDE, cM + cN);
    VecSetUp(rhs);

    Vec x;

    VecDuplicate(rhs, &x);

    double * arrayg1;
    double * arrayg2;
    double * arrayrhs;

    VecGetArray(g1, &arrayg1);
    VecGetArray(g2, &arrayg2);
    VecGetArray(rhs, &arrayrhs);

    for (int i=0; i<cM; i++){
        arrayrhs[i] = arrayg1[i]; 
    }

    for (int i=0; i<cN; i++){
        arrayrhs[i+cM] = arrayg2[i]; 
    }

    VecRestoreArray(g1, &arrayg1); 
    VecRestoreArray(g2, &arrayg2); 
    VecRestoreArray(rhs, &arrayrhs);

    KSPSolve(ksp, rhs, x);

// Check computed system

    //cout << "Boundary dof size : " << bndryDarcy.size() << endl;

   //for (const auto& it: bndryDarcy){
   //    cout << it.first << endl;
   //}

    const char *checkA = "MatrixCheckA.dat";

    // Write A matrix
    WriteMat(reducedsys->M,checkA);

    const char *checkB = "MatrixCheckB.dat";

    // Write B matrix
    WriteMat(reducedsys->B,checkB);

    // Write right two right hand side vectors
    const char *checkg1 = "MatrixCheckg1.dat";

    WriteVec(g1, checkg1);   

    WriteMat(reducedsys->B,checkB);

    const char *checkg2 = "MatrixCheckg2.dat";

    WriteVec(g2, checkg2);   

    //MatView(matrix->Bd, PETSC_VIEWER_STDOUT_WORLD);
    //MatView(reducedsys->B, PETSC_VIEWER_STDOUT_WORLD);
    //MatView(reducedsys->M, PETSC_VIEWER_STDOUT_WORLD);
    //VecView(g1, PETSC_VIEWER_STDOUT_WORLD);
    //VecView(g2, PETSC_VIEWER_STDOUT_WORLD);
    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    cout << "here !" << endl;

    // Test inexect Uzawa iteration algorithm
    linearSys * ls = (linearSys *)malloc(sizeof(linearSys));

    PetscCall(MatConvert(BT, MATSAME, MAT_INITIAL_MATRIX, &ls->B));
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

    PreconditionedUzawa(ls, 10e-10, 5);

    VecView(ls->x, PETSC_VIEWER_STDOUT_WORLD);
    VecView(ls->y, PETSC_VIEWER_STDOUT_WORLD);

// ====================================================================================================================================
    // Clear used objects
    DMDAVecRestoreArray(dmu,localu,&lu);
    DMRestoreLocalVector(dmu, &localu); 

    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    //VecDestroy(&source);
    DMDestroy(&dm);
    DMDestroy(&dmu);

    PetscFinalize();

    return 0;
}
