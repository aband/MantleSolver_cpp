#include <iostream>
#include <petsc.h>
#include "integral.h"

#include "stencil.h"
#include "util.h"
#include "polynomial.h"
#include "input.h"
#include "reconstruction.h"
//#include <adolc/adolc.h>

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

double func(vertex& point, const vector<double>& param){
	 //if (point[0]<-1.0/param[0]){
//		  return point[0]*point[0]+point[1]*point[1];
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1);
//	 } else {
//		  return point[0]*point[0]*point[1]*point[1] + 1.0;
//	     return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1) + 1;
//	 }

    return sin(point[0]*3.0)+cos(point[1]/2.0) + point[0]*(point[1]+1);
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

    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // ==========================================================================================================================

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 3, N = 3;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    // Create data management object
    DM    dm;
    Vec   fullmesh;
    const int stencilWidth = 5;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidth, NULL, NULL, &dm);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dm);               CHKERRQ(ierr);
    ierr = DMSetUp(dm);                        CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm, &fullmesh);CHKERRQ(ierr); 

    double L = 2.0, H = 2.0;
    double xstart = -1.0, ystart = -1.0;
    ierr = PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL); CHKERRQ(ierr);

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

    cout << "Mesh Created. To check full mesh, rerun with -printmesh 1 " << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // ==========================================================================================================================

    // Contain defined mesh in vector container.
    // and verify it.
    vector< valarray<double> > mesh;
    
    ReadMeshPortion(dm, &fullmesh, mesh);

    //cout << "Converted c array of local mesh into vector container c++ " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // ====================================================================================================================================

    DM dmu;

    int cell_ghost = 3;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, cell_ghost, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu);               CHKERRQ(ierr);
    ierr = DMSetUp(dmu);                        CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    // Initialize with oblique data for Burgers equation 
    //ObliqueBurgers(dm,dmu,&fullmesh,&globalu,Initial_Condition);
    SimpleInitialValue(dm,dmu,&fullmesh,&globalu,func);

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
    AssignValuesMeshInfo(mi,dm,dmu); 

    // Assign local mesh and local values to mi
    mi.lmesh = mesh;
    mi.localVals = lu;

    indice start = {M/2,N/2};
    vector<indice> targetCell;
    targetCell.push_back({0,0});
    targetCell.push_back({1,0}); 
    vertex center = {0.0,0.0};

    MLWENO::stencil<indice> third(3,3);

    for (int j=0; j<3; j++){
    for (int i=0; i<3; i++){
        third(i,j) = {i-1,j-1};
    }}

    MLWENO::stencil<indice> second(2,2);

    for (int j=0; j<2; j++){
    for (int i=0; i<2; i++){
        second(i,j) = {i,j};
    }}

    //MLWENO::stencilPolynomial* sp = new MLWENO::stencilPolynomial(start,center,targetCell);

    //sp->SetStencilPolynomials(mi, second);

    //sp->printCoef(); 
    
    //cout << sp->eval(0.0,0.0) << " " <<func(center,{1,1}) << endl;

    // Test class of reconstruction
    MLWENO::reconstruction * rptr = new MLWENO::reconstruction();

    int stencil3[2] = {3,3};
    int stencil2[2] = {2,2};
    vector<int*> stencilSizes = {stencil3,stencil2,stencil2,stencil2,stencil2};
    vector<indice> shifts = {{-1,-1},{-1,-1},{0,0},{-1,0},{0,-1}};

    rptr->AddStencil(stencilSizes, shifts);

    rptr->CreateStencilPolynomials(start, center, targetCell, mi);

    rptr->Update(mi);

    //rptr->PrintNonLinWgts();

    //rptr->PrintSmoothnessIndic();

    //rptr->PrintStencils();

    //cout << rptr->Eval(0.0,0.0) << " " << func(center,{(double)M,(double)N}) << endl;

    rptr->Clear();

    // Test single level reconstruction

    //MLWENO::singleLevelReconstruction * slrPtr = new MLWENO::singleLevelReconstruction(3,3);

    //slrPtr->CreateStencilPolynomials(mi);

    //indice test = {0,0};
    //slrPtr->CheckStencilPolynomials(mi,test);

    // Test multi level reconstruction
    MLWENO::multiLevelReconstruction * mlrPtr = new MLWENO::multiLevelReconstruction(mi,2,2);
    mlrPtr->AddLevel(mi,3,3);
    mlrPtr->AddLevel(mi,2,3);
    mlrPtr->AddLevel(mi,3,2);

    vector<vector<indice>> reconstMethod {{{-1,0},{-1,-1,},{0,-1},{0,0}},{{-1,-1}},{{-1,-1},{0,-1}},{{-1,-1},{-1,0}}};

    //vector<vector<indice>> reconstMethod {{{-1,0},{-1,-1,},{0,-1},{0,0}},{{-1,-1}}};

    mlrPtr->AddReconstMethod(reconstMethod);

    //mlrPtr->GetInfo();

    mlrPtr->UpdateTwoStageNonLinearWgts(mi);

    mlrPtr->PrintNonLinearWgts(mi);

    //mlrPtr->PrintSmoothnessIndicator(mi);

    // ====================================================================================================================================
    // Clear used objects
    DMDAVecRestoreArray(dmu,localu,&lu);
    DMRestoreLocalVector(dmu, &localu); 

    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    DMDestroy(&dm);
    DMDestroy(&dmu);

    return 0;
}
