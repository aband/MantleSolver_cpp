#include "driver.h"

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

    // Infinitly smooth test case
    //return sin(point[0] + point[1] + 0.1);

    // Quadratic test case
    //return point[0]*point[0] + point[1]*point[1];

    // Linear test case
    //return point[0] + point[1];

    // Constant test case
    //return 0.5;
}

int main(int argc, char **argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Initialize mesh ==========================================================
    int M=3, N=3;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL));

    DM dmMesh;

    int stencilWidthMesh = 5;

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, &dmMesh));
    PetscCall(DMSetFromOptions(dmMesh));              
    PetscCall(DMSetUp(dmMesh));

    Vec globalmesh;

    PetscCall(DMCreateGlobalVector(dmMesh, &globalmesh));

    // Physical domain
    double L = 2.0, H = 2.0;
    double xstart = -1.0, ystart = -1.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    // If running test on the single stencil for convergence study
    int singleStencilTest = 0;
    double scale = 1;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-single", &singleStencilTest, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL, "-scale", &scale, NULL));

    if (singleStencilTest){
        L = L/scale;
        H = H/scale;
        xstart = -L/2.0;
        ystart = -H/2.0;
    }

    // Define Intermediate data structure
    MeshParam mp;
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    int meshtype = 0; 

    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshtype,NULL));
    switch(meshtype){
        case 0: CreateFullMesh(dmMesh, &globalmesh, &mp); break;
        case 1: LogicRectMesh(dmMesh, &globalmesh, &mp);  break;
        case 2: RefineMesh(dmMesh, &globalmesh, &mp);
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    // Print mesh or not
    int printmesh=0;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-printmesh",&printmesh,NULL));
    if(printmesh){ 
        VecView(globalmesh, PETSC_VIEWER_STDOUT_WORLD);
        PrintFullMesh(dmMesh, &globalmesh);
    }

    // Setting initial solution
	 DM dmu; 
    int stencilWidthU = 3;

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 
    stencilWidthU, NULL, NULL, &dmu));
    PetscCall(DMSetFromOptions(dmu));              
    PetscCall(DMSetUp(dmu));     

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));


    // Set up global vector holding cell averaged solution
    Vec globalu;
    PetscCall(DMCreateGlobalVector(dmu,&globalu));

    // Set up initial values for cell averaged solution
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalu, func); 

    // Scattering global solution into local pieces
    Vec localu; 
    PetscCall(DMGetLocalVector(dmu, &localu));

    PetscCall(DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu));
    PetscCall(DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu));


    // Set up MeshInfo struct ==================================================
    MeshInfo mi;

    DMDAVecGetArray(dmu, localu, &mi.localVals);

    ReadMeshPortion(dmMesh, &globalmesh, mi.lmesh);

    // Assign meshinfo after mesh added to meshinfo
    AssignValuesMeshInfo(mi,dmMesh,dmu); 

    Driver * driPtr = new Driver(&mi);


    // Finialie the program ====================================================
    PetscCall(DMDAVecRestoreArray(dmu,localu,&mi.localVals));
    PetscCall(DMRestoreLocalVector(dmu, &localu));

    PetscCall(VecDestroy(&globalu));
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    PetscFinalize();

    return 0;
}
