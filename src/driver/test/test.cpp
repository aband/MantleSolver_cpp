#include "driver.h"
#include "mlwenouse.h"

double tmpfunc(const vertex& point, const vector<double>& param){
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

bool left_boundary(const indice& globalCell,
                   const MeshInfo& mi){

    if (globalCell[0] == 0 &&
        globalCell[1] > 0 && 
        globalCell[1] < mi.MPIglobalCellSize[1]-1) {
        return true;
    } else {
        return false;
    }
}

bool interior(const indice& globalCell, 
              const MeshInfo& mi){
   if (globalCell[0] > 0 && globalCell[0] < mi.MPIglobalCellSize[0]-1 &&
       globalCell[1] > 0 && globalCell[1] < mi.MPIglobalCellSize[1]-1){
       return true;
   } else {
       return false;
   }
}

bool edge(const indice& globalCell,
          const MeshInfo& mi){

   // return four edges

   if (// left edge
       (globalCell[0] == 0 && 
        globalCell[1] != 0 && 
        globalCell[1] != mi.MPIglobalCellSize[1]-1) ||
       // right edge
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && 
        globalCell[1] != 0 && 
        globalCell[1] != mi.MPIglobalCellSize[1]-1) ||
       // bottom edge
       (globalCell[1] == 0 && 
        globalCell[0] != 0 && 
        globalCell[0] != mi.MPIglobalCellSize[0]-1) ||
       // top edge
       (globalCell[1] == mi.MPIglobalCellSize[1]-1 && 
        globalCell[0] != 0 && 
        globalCell[0] != mi.MPIglobalCellSize[0]-1) 
      ){
       return true;
   } else {
       return false;
   }

}

bool corner(const indice& globalCell, 
            const MeshInfo& mi){

   // return four corners

   if ((globalCell[0] == 0 && globalCell[1] == 0) ||
       (globalCell[0] == 0 && globalCell[1] == mi.MPIglobalCellSize[1]-1) ||
       (globalCell[0] == 0 && globalCell[1] == 0) ||
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && globalCell[1] == mi.MPIglobalCellSize[1]-1) 
      ){
       return true;
   } else {
       return false;
   }

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
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalu, {0.0,0.0},tmpfunc); 

    // Scattering global solution into local pieces
    Vec localu; 
    PetscCall(DMGetLocalVector(dmu, &localu));

    PetscCall(DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu));
    PetscCall(DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu));

    // Set up Driver class =====================================================
    Driver * drivPtr = new Driver();

    DMDAVecGetArray(dmu, localu, &drivPtr->mi.localVals);

    ReadMeshPortion(dmMesh, &globalmesh, drivPtr->mi.lmesh);

    // Assign meshinfo after mesh added to meshinfo
    AssignValuesMeshInfo(drivPtr->mi,dmMesh,dmu); 

    drivPtr->UseWeno();
    drivPtr->AddLevel(1,1);
    drivPtr->AddLevel(2,2);
    drivPtr->AddLevel(3,3);
    drivPtr->AddLevel(4,4);
    drivPtr->AddLevel(5,5);

    // =========================================================================
    MLWENO::MLWENOPrepare * mlpPtr = new MLWENO::MLWENOPrepare();

    mlpPtr->AddLevel(drivPtr->mi,1,1);
    mlpPtr->AddLevel(drivPtr->mi,2,2);
    mlpPtr->AddLevel(drivPtr->mi,3,3);
    mlpPtr->AddLevel(drivPtr->mi,4,4);
    mlpPtr->AddLevel(drivPtr->mi,5,5);

    mlpPtr->UpdateSmoothnessIndic(drivPtr->mi);

    MLWENO::MLWENOUse * mluse = new MLWENO::MLWENOUse(); 

    // Add reconstruction levels to interior cells
    mluse->AddMLWENOLevel("interior",{"(3,3)","(5,5)"}, mlpPtr);

    mluse->AssignWENOStencils("interior","(3,3)",{{-2,0},{-2,-2},{0,0},{0,-2}});
    mluse->AssignWENOStencils(0,"(5,5)",{{-2,-2}});

    mluse->UpdateNonLinearWgts(drivPtr->mi, "interior", "two_stage", interior);

    // Add reconstruction levels to boundary cells
    mluse->AddMLWENOLevel("left_boundary",{"(3,3)","(1,1)"}, mlpPtr);

    mluse->AssignWENOStencils("left_boundary","(3,3)",{{0,0},{0,-2}});
    mluse->AssignWENOStencils(1,"(1,1)",{{0,0}});

    mluse->UpdateNonLinearWgts(drivPtr->mi, "left_boundary", "one_stage", left_boundary);

    mluse->PrintNonLinearWgts("interior",drivPtr->mi);

    /**!
     * L_1 error
     */
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;
    double work = 0.0;

    indice cell {M/2,N/2}; 

    vector<vertex> corner = extractCorners(drivPtr->mi, cell); 

    PetscCall(PetscPrintf(PETSC_COMM_SELF,"The point-wise error at point (%f, %f) is %.3e \n",
                                           0.0,0.0,abs(mluse->Evaluate({0,0},{M/2,N/2},drivPtr->mi, "interior")-tmpfunc({0,0},{0,0}))));

    for (int g=0; g<gpf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g], corner);
        double jac = abs(GaussJacobian(gpf[g], corner));
        work += gwf[g]*abs(tmpfunc(mapped,{0,0}) - mluse->Evaluate(mapped, cell, drivPtr->mi, "interior"))*jac;
    }

    PetscCall(PetscPrintf(PETSC_COMM_SELF,"The L1 Error at cell (%d, %d) is %.3e \n",
                                           cell[0],cell[1],work));

    // Finialize the program ====================================================
    PetscCall(DMDAVecRestoreArray(dmu,localu,&drivPtr->mi.localVals));
    PetscCall(DMRestoreLocalVector(dmu, &localu));

    PetscCall(VecDestroy(&globalu));
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    PetscFinalize();

    return 0;
}
