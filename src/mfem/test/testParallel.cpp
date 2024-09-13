#include <iostream>
#include <fstream>
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
#include "passemble.h"
#include "preconst.h"
#include "psolve.h"
#include <ctime>
#include <chrono>

extern "C"{
#include "mesh.h"
#include "output.h"
#include "cgns_io.h"
}

using namespace std;

int main(int argc, char ** argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Parallel test \n"));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));

    // =================================================================================

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 2, N = 2;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    // Create data management object
    DM    dm;
    Vec   fullmesh;

    const int stencilWidth = 2;

    int partition = 0;
    PetscOptionsGetInt(NULL, NULL, "-part", &partition, NULL);

    int m = 1,n = 1;
    PetscOptionsGetInt(NULL, NULL, "-mym", &m, NULL);
    PetscOptionsGetInt(NULL, NULL, "-myn", &n, NULL);

    PetscInt *lx;
    PetscInt *ly;

    PetscMalloc1(m, &lx);
    PetscMalloc1(n, &ly);

    if (partition == 1){
    if (m*n == 0 || m*n != size){
        PetscPrintf(PETSC_COMM_WORLD,"Please provide valid partitions\n");
        return 0;
    }}

    // Careate partition array
    int mpr = M/m;
    int npr = N/n;

    for (int i=0; i<m-1; i++){
        lx[i] = mpr;
    }

    for (int i=0; i<n-1; i++){
        ly[i] = npr;
    }

    lx[m-1] = M - mpr*(m-1);
    ly[n-1] = N - npr*(n-1);

    if (partition == 0){
        ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidth, NULL, NULL, &dm);CHKERRQ(ierr);
        ierr = DMSetFromOptions(dm);               CHKERRQ(ierr);
        ierr = DMSetUp(dm);                        CHKERRQ(ierr);
    } else {
        ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, m, n, 2, stencilWidth, lx, ly, &dm);CHKERRQ(ierr);
        ierr = DMSetFromOptions(dm);               CHKERRQ(ierr);
        ierr = DMSetUp(dm);                        CHKERRQ(ierr);
    }

    ierr = DMCreateGlobalVector(dm, &fullmesh);CHKERRQ(ierr); 

    PhysProperty * physproperty = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(physproperty);

    double physscale = physproperty->L0/physproperty->l0;
    double L = 2*physscale, H = 1*physscale;
    double xstart = -1*physscale, ystart = -1.0001*physscale;
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

    int cell_ghost = 1;

  
    if (partition == 0){ 
        PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
                               DM_BOUNDARY_PERIODIC, 
                               DM_BOUNDARY_PERIODIC, 
                               DMDA_STENCIL_BOX, 
                               M,N, 
                               PETSC_DECIDE, 
                               PETSC_DECIDE, 
                               1, cell_ghost, NULL, NULL, &dmu));
        PetscCall(DMSetFromOptions(dmu)); 
        PetscCall(DMSetUp(dmu)); 

    } else {

        PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
                               DM_BOUNDARY_PERIODIC, 
                               DM_BOUNDARY_PERIODIC, 
                               DMDA_STENCIL_BOX, 
                               M,N, 
                               m,n, 
                               1, cell_ghost, lx, ly, &dmu));
        PetscCall(DMSetFromOptions(dmu)); 
        PetscCall(DMSetUp(dmu)); 

        PetscCall(DMView(dmu, PETSC_VIEWER_STDOUT_WORLD));
    }

    // === Evenly divided physical domain

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

    //ParallelAssembleTest();

    // =================================================================================
    // Declare classes for shape functions
    basis * basis_   = new basis();
    Hdivmixed * hdiv = new Hdivmixed();
    BRMixed * br     = new BRMixed();

    br->ComputeTotalDOF(mi);
    hdiv->ComputeTotalDOF(mi);

    // Mark boundary values
    bndryVal bndryStokesEssen;
    bndryVal bndryStokesNatur;
    bndryVal bndryDarcyEssen;
    bndryVal bndryDarcyNatur;

    MarkBndryDOFStokes(bndryStokesEssen, bndryStokesNatur, mi, *basis_, *br, physproperty);
    MarkBndryDOFDarcy(bndryDarcyEssen, bndryDarcyNatur, mi, *basis_, *hdiv, physproperty);

    // Create linear system
    ReducedSys * reducedDarcy = (ReducedSys *)malloc(sizeof(ReducedSys));
    ReducedSys * reducedStokes = (ReducedSys *)malloc(sizeof(ReducedSys));

    Mat K;

    int bndryDOFStokes = 0.0;
    int bndryDOFDarcy  = 0.0;

    int * refArrayStokes = new int[br->getDOF()];
    int * refArrayDarcy  = new int[hdiv->getDOF()];

    CreateRefMap(*br, refArrayStokes, mi, &bndryDOFStokes);
    CreateRefMap(*hdiv, refArrayDarcy, mi, &bndryDOFDarcy);

    ParallelMatrixAssemble(mi, *basis_, physproperty, bndryStokesEssen, reducedStokes, 
                                                      bndryDarcyEssen,  reducedDarcy, 
                           &K, *br, *hdiv , refArrayStokes, refArrayDarcy, bndryDOFStokes, bndryDOFDarcy);

    int nelem = M*N;
    CreateLinearSys(reducedStokes, nelem);
    CreateLinearSys(reducedDarcy , nelem); 

    // ======================================================================
    // Create Reduced system using previous routine
/*
    // Assemble matrices, create full system
    System * system = (System *)malloc(sizeof(System));

    SerialMatrixAssembleBlock(mi, *basis_, *hdiv, *br, physproperty, system);

    ReducedSys * reducedDarcySerial = (ReducedSys *)malloc(sizeof(ReducedSys));
    ReducedSys * reducedStokesSerial = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateReducedSerial(reducedDarcySerial, &system->Ad, &system->Bd, &system->sourceDarcy, bndryDarcyEssen);
    CreateReducedSerial(reducedStokesSerial, &system->As, &system->Bs, &system->sourceStokes, bndryStokesEssen);

   CreateLinearSys(reducedStokesSerial, nelem);
   CreateLinearSys(reducedDarcySerial, nelem);

   MatAXPY(reducedDarcySerial->M, -1.0,reducedDarcy->M,SUBSET_NONZERO_PATTERN);  
   MatAXPY(reducedStokesSerial->M, -1.0,reducedStokes->M,SUBSET_NONZERO_PATTERN); 
   MatView(reducedDarcySerial->M, PETSC_VIEWER_STDOUT_WORLD);
   MatView(reducedStokesSerial->M, PETSC_VIEWER_STDOUT_WORLD);

   MatAXPY(reducedDarcySerial->B, -1.0,reducedDarcy->B,SUBSET_NONZERO_PATTERN);  
   MatAXPY(reducedStokesSerial->B, -1.0,reducedStokes->B,SUBSET_NONZERO_PATTERN); 
   MatView(reducedDarcySerial->B, PETSC_VIEWER_STDOUT_WORLD);
   MatView(reducedStokesSerial->B, PETSC_VIEWER_STDOUT_WORLD);

   MatAXPY(system->Cd, -1.0,reducedDarcy->C,SUBSET_NONZERO_PATTERN);  
   MatAXPY(system->Cs, -1.0,reducedStokes->C,SUBSET_NONZERO_PATTERN); 
   MatView(system->Cd, PETSC_VIEWER_STDOUT_WORLD);
   MatView(system->Cs, PETSC_VIEWER_STDOUT_WORLD);

   VecAXPY(reducedStokesSerial->F, -1.0, reducedStokes->F);
   VecView(reducedStokesSerial->F, PETSC_VIEWER_STDOUT_WORLD);

   VecAXPY(reducedDarcySerial->F, -1.0, reducedDarcy->F);
   VecView(reducedDarcySerial->F, PETSC_VIEWER_STDOUT_WORLD);
*/
//    MatConvert(reducedStokesSerial->Kg, MATAIJ, MAT_INITIAL_MATRIX, &reducedStokes->Kg);
//    MatConvert(reducedStokesSerial->M, MATAIJ, MAT_INITIAL_MATRIX, &reducedStokes->M);

    // ======================================================================

    // Create the coupled saddle point system
    // Attention, stokes system should be used as the first input
    ReducedSys * Result = (ReducedSys *)malloc(sizeof(ReducedSys));

    CreateCoupledSystem(reducedStokes, reducedDarcy, Result, &K);

    // Solve with Uzawa solver
    int maxIter = 1;
    PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL);
    double tolUzawa = 10e-7;
    PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL);

    auto start = std::chrono::system_clock::now();
    CoupledUzawa(Result, tolUzawa, maxIter);    
    auto end = std::chrono::system_clock::now();
 
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    PetscCall(PetscPrintf(PETSC_COMM_SELF, "Current rank is %d : Computation time is : %f s \n", rank, elapsed_seconds.count()));

    // Create flow velocity at cell centroid
    int nelemloc = mi.MPIlocalCellSize[0]*mi.MPIlocalCellSize[1];
    double * ux = (double *)malloc(sizeof(double)*nelemloc);
    double * uy = (double *)malloc(sizeof(double)*nelemloc);

    double * vx = (double *)malloc(sizeof(double)*nelemloc);
    double * vy = (double *)malloc(sizeof(double)*nelemloc);

    Vec stokesv;
    Vec darcyv;

    PetscCall(VecNestGetSubVec(Result->x, 0, &stokesv));
    PetscCall(VecNestGetSubVec(Result->x, 1, &darcyv));

    Vec destStokes_sol, destStokes_g;
    Vec destDarcy_sol, destDarcy_g;

    SolScatAll(&stokesv, &reducedStokes->g, 
               &destStokes_sol, &destStokes_g);  

    SolScatAll(&darcyv, &reducedDarcy->g, 
               &destDarcy_sol, &destDarcy_g);  

    CGNSPrepareParallel(&destStokes_sol, &destStokes_g, refArrayStokes, mi,
                        ux, uy, *br, *basis_);

    CGNSPrepareParallel(&destDarcy_sol, &destDarcy_g, refArrayDarcy, mi,
                        vx, vy, *hdiv, *basis_);

    // CGNS output of hdf5 file
    char stokesfile[] = "stokes.cgns";   
    CgnsArrayOutput(dm,&fullmesh,ux,uy,mi.MPIlocalCellStart[0],
                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
                    mi.MPIlocalCellSize[1],stokesfile);

    char darcyfile[] = "darcy.cgns";    	
    CgnsArrayOutput(dm,&fullmesh,vx,vy,mi.MPIlocalCellStart[0],
                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
                    mi.MPIlocalCellSize[1],darcyfile);

    // Finalize Petsc code
    VecDestroy(&destStokes_sol);
    VecDestroy(&destStokes_g);
    VecDestroy(&destDarcy_sol);
    VecDestroy(&destDarcy_g);

    DMDAVecRestoreArray(dmu,localu,&lu);
    DMRestoreLocalVector(dmu, &localu); 

    VecDestroy(&fullmesh);
    VecDestroy(&globalu);

    DMDestroy(&dm);
    DMDestroy(&dmu);

    free(refArrayStokes);
    free(refArrayDarcy);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Run in parallel successfully! ^_^\n"));

    PetscFinalize();

    return 0;
}
