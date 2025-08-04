#include "petsc.h"
#include "input.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

#include "tensorstencilpoly.h"
#include <chrono>

#include "error.h"

double func(const vertex& point,
            const vector<double>& param){

    if (point[0] < param[0]) {

    return sin(point[0])*cos(point[1]);

    } else {

    return sin(point[0])*cos(point[1]) + 0.1;

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

    // Pick a stencil
    int startM = M/2-2;
    int startN = N/2-2;

    // Define tensor product stencil polynomial
    tensorstencilpoly stenpoly = tensorstencilpoly(2);

    stenpoly.setCoef(mi, startM, startN);
    stenpoly.printCoef();

    cout << endl;
    for (int cell=0; cell<9; cell++){
        cout << stenpoly.eval(0.46,0.4505,cell) << "  ";
    }
    cout << endl;

    double val[9] = {0};

    //for (int ncell = 0; ncell<9; ncell ++){
    int ncell = 1;
    stenpoly.der(2,2,ncell,val, 0.46, 0.4505); 

    cout << endl;
    for (int j=0; j<3; j++){
        for (int i=0; i<3; i++){
            cout << val[j*3+i] << "   " ;
        } cout << endl;
    }

    cout << endl;
    //}

    stenpoly.setSigma();
    stenpoly.printSigmaBase();

    // Set full stencils
    cout << M << "  " << N << endl;
    vector<tensorstencilpoly *> allsten;
    allsten.resize(M*N);

	 auto start = std::chrono::steady_clock::now();
    for (int j=0; j<N-4; j++){
    for (int i=0; i<M-4; i++){
			int s = j*M+i;
        allsten.at(s) = new tensorstencilpoly(5);
        allsten.at(s)->setCoef(mi,i,j);
		  allsten.at(s)->setSigma();
    }}
	 auto end = std::chrono::steady_clock::now();
	 auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    cout << "Time Test: " << duration.count() << " ms." << endl;

    // Reconstruction test
    Vec globalvec, localvec;
    double ** locvals;

    PetscCall(DMCreateGlobalVector(dmu, &globalvec));

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalvec, {0.0}, func);
//VecView(globalvec, PETSC_VIEWER_STDOUT_WORLD);
    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvec)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalvec, INSERT_VALUES, localvec));
    PetscCall(DMGlobalToLocalEnd(dmu, globalvec, INSERT_VALUES, localvec));

    PetscCall(DMDAVecGetArray(dmu, localvec, &locvals));

    printexactsol(mi, 0, func, 1, true, {0.0});

    // =================================================================
    DMDAVecRestoreArray(dmu,localvec,&locvals);
    DMRestoreLocalVector(dmu, &localvec); 

    PetscCall(VecDestroy(&globalvec));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));


    return 1;
}
