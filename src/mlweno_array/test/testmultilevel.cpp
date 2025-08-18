#include "petsc.h"
#include "input.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

#include "tensorstencilpoly.h"
//#include "reconstruction.h"
#include "multilevel.h"
#include <chrono>

#include "error.h"

double func(const vertex& point,
            const vector<double>& param){

    if (point[0] < 0.5) {

    return sin(point[0])*cos(point[1]);

    } else {

    return sin(point[0])*cos(point[1]) ;

    }
}

int main(int argc, char ** argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    int M = 30, N = 30;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    double L = 1, H = 1;
    //double xstart = -L/2, ystart = -H/2;
    double xstart = 0.0, ystart = 0.0;
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

    double dscale = 1.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-scale",&dscale,NULL));

    L/=dscale;
    H/=dscale;

    double vscale = 1.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-vscale",&vscale,NULL));

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

    mi.L = L;
    mi.H = H;

/*
    vector<tensorstencilpoly> sten3;
    sten3.resize((M-2)*(N-2));
    for (int j=0; j<N-2; j++){
    for (int i=0; i<M-2; i++){
        int s = j*(M-2)+i;
        sten3.at(s) = tensorstencilpoly(2,3,3);
        sten3.at(s).setCoef(mi,i,j);
		  sten3.at(s).setSigma();
		  sten3.at(s).startx = i;
		  sten3.at(s).starty = j;
    }}

    sten3.at(0).printCoef();
*/

    // =================================================================
    multilevel ml = multilevel();
    ml.addLevel(1,"(2,2)",2,2,mi);
    ml.addLevel(2,"(3,3)",3,3,mi);
    ml.addLevel(3,"(4,4)",4,4,mi);
    ml.addLevel(4,"(5,5)",5,5,mi);

    cout << "Levels created ..." << endl;

    ml.printInfo();

    unordered_map<std::string, vector<indice>> my53method;
    my53method.insert(std::make_pair<std::string, vector<indice>>("(5,5)", {{-2,-2}}));
    my53method.insert(std::make_pair<std::string, vector<indice>>("(3,3)", {{-2,-2},{-2,0},{0,0},{0,-2}}));

    unordered_map<std::string, indice> mysize;
	 mysize.insert(std::make_pair<std::string, indice>("(5,5)", {5,5}));
	 mysize.insert(std::make_pair<std::string, indice>("(3,3)", {3,3}));

    ml.my_recon.resize(M*N);

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
      	int s = j*M+i;
        ml.my_recon.at(s) = mlreconstruction();
        ml.my_recon.at(s).prepare(my53method, mysize, {i,j}, true, mi);
		  cout <<endl << "Cell : "<< i << ", "<< j << endl;
		  ml.my_recon.at(s).printInfo();
    }}

    // Reconstruction test
    Vec globalvec, localvec;
    double ** locvals;

    PetscCall(DMCreateGlobalVector(dmu, &globalvec));

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalvec, {0.0}, func);
    // VecView(globalvec, PETSC_VIEWER_STDOUT_WORLD);
    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvec)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalvec, INSERT_VALUES, localvec));
    PetscCall(DMGlobalToLocalEnd(dmu, globalvec, INSERT_VALUES, localvec));

    PetscCall(DMDAVecGetArray(dmu, localvec, &locvals));



    // =================================================================
    DMDAVecRestoreArray(dmu,localvec,&locvals);
    DMRestoreLocalVector(dmu, &localvec); 

    PetscCall(VecDestroy(&globalvec));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    return 1;
}
