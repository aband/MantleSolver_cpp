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

    Vec source;

    br->ComputeTotalDOF(mi);

    hdiv->ComputeTotalDOF(mi);

    DMCreateGlobalVector(dm, &source);

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

    CreateReducedSerial(reducedsys, &matrix->Ad, bndryDarcy);

    // Check mat size
    int cM, cN;
    PetscCall(MatGetSize(reducedsys->Kg, &cM, &cN));
    cout << cM << " " << cN << endl;

    Vec g;
    PetscCall(VecDuplicate(reducedsys->g, &g));
    PetscCall(VecGetSize(g, &cN));
    cout << cN << endl;

    PetscCall(MatMult(reducedsys->Kg, reducedsys->g, g));


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
