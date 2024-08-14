#include <petsc.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include "integral.h"
#include "phase.h"
#include "param.h"
#include "input.h"
#include "util.h"

// MFEM parameter header file
#include "myFunc.h"
#include "passemble.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "bndry.h"
#include "preconst.h"
#include "psolve.h"

// MLWENO parameter header file
#include "driver.h"
#include "mlwenouse.h"

extern "C"{
#include "mesh.h"
#include "output.h"
//#include "cgns_io.h"
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

    // Create phase class containing constant physic attributes and 
    // phase behavior package.
    Phase * phase = new Phase();

    phase->pp = (PhysProperty *)malloc(sizeof(PhysProperty));

    AssignPhyProperties(phase->pp);

    phase->pPtr = new phaseState();

    // Physical domain
    double physscale = phase->pp->L0/phase->pp->l0;
    double L = 2*physscale, H = 1*physscale;
    double xstart = -1*physscale, ystart = -1.0001*physscale;

    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    DM dmMesh;

    int stencilWidthMesh = 5;

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, &dmMesh));
    PetscCall(DMSetFromOptions(dmMesh));              
    PetscCall(DMSetUp(dmMesh));

    Vec globalmesh;

    PetscCall(DMCreateGlobalVector(dmMesh, &globalmesh));

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
    // ==========================================================================

    // Create initial (C,H) distribution pair
    // Calculate volume fraction of fluid (porosity)
    PorosityOut(xstart, ystart, L, H, 20, phase);      

    // Solve for velocity with finite element solver
    vector<valarray<double>> mesh;
    ReadMeshPortion(dmMesh, &globalmesh, mesh);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    Vec localu; 
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu);

    // It can be changed later to not be double
    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    MeshInfo mi;
    mi.lmesh     = mesh;
    mi.localVals = lu;

    AssignValuesMeshInfo(mi, dmMesh, dmu);

    // Get shape functions;
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

    MarkBndryDOFStokes(bndryStokesEssen, bndryStokesNatur, mi, *basis_, *br, phase->pp);
    MarkBndryDOFDarcy(bndryDarcyEssen, bndryDarcyNatur, mi, *basis_, *hdiv, phase->pp);

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

    ParallelMatrixAssemble(mi, *basis_, phase, bndryStokesEssen, reducedStokes, 
                                               bndryDarcyEssen,  reducedDarcy, 
                           &K, *br, *hdiv , refArrayStokes, refArrayDarcy, bndryDOFStokes, bndryDOFDarcy);

    int nelem = M*N;
    CreateLinearSys(reducedStokes, nelem);
    CreateLinearSys(reducedDarcy , nelem); 

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

    // Sequential visual output
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
//    char stokesfile[] = "stokes.cgns";   
//    CgnsArrayOutput(dmMesh,&globalmesh,ux,uy,mi.MPIlocalCellStart[0],
//                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
//                    mi.MPIlocalCellSize[1],stokesfile);

//    char darcyfile[] = "darcy.cgns";    	
//    CgnsArrayOutput(dmMesh,&globalmesh,vx,vy,mi.MPIlocalCellStart[0],
//                    mi.MPIlocalCellSize[0], mi.MPIlocalCellStart[1],
//                    mi.MPIlocalCellSize[1],darcyfile);

    // Transport ================================================================
    // Create levels for ml-weno 
    Driver * drivPtr = new Driver();
    DMDAVecGetArray(dmu, localu, &drivPtr->mi.localVals);
    ReadMeshPortion(dmMesh, &globalmesh, drivPtr->mi.lmesh);

    // Assign mesh information after mesh added to meshInfo
    AssignValuesMeshInfo(drivPtr->mi, dmMesh, dmu);

    // We have five different levels in ml-weno 
    drivPtr->UseWeno();
	 drivPtr->AddLevel(1,1);
	 drivPtr->AddLevel(2,2);
    drivPtr->AddLevel(3,3);
    drivPtr->AddLevel(4,4);
    drivPtr->AddLevel(5,5);

    MLWENO::MLWENOPrepare * mlpPtr = new MLWENO::MLWENOPrepare();

    mlpPtr->AddLevel(drivPtr->mi,1,1);
    mlpPtr->AddLevel(drivPtr->mi,2,2);
    mlpPtr->AddLevel(drivPtr->mi,3,3);
    mlpPtr->AddLevel(drivPtr->mi,4,4);
    mlpPtr->AddLevel(drivPtr->mi,5,5);

    mlpPtr->UpdateSmoothnessIndic(drivPtr->mi);

    // Two instance of mlweno usage, advection and diffusion
    // Advection mlweno use (3,3) and (2,2) reconstruction
    MLWENO::MLWENOUse * mluseAdv = new MLWENO::MLWENOUse(); 

    mluseAdv->AddMLWENOLevel("interior",{"(3,3)","(2,2)"}, mlpPtr);

    mluseAdv->AssignWENOStencils(0,"(2,2)",{{-1,0},{-1,-1},{0,0},{0,-1}});
    mluseAdv->AssignWENOStencils(0,"(3,3)",{{-1,-1}});

    mluseAdv->UpdateNonLinearWgts(drivPtr->mi, "interior", "two_stage", interior);

    // Diffusion mlweno use (5,5) and (3,3) reconstruction
    MLWENO::MLWENOUse * mluseDif = new MLWENO::MLWENOUse();

    mluseDif->AddMLWENOLevel("interior",{"(5,5)","(3,3)"}, mlpPtr);

    mluseDif->AssignWENOStencils(0,"(3,3)",{{-2,0},{-2,-2},{0,0},{0,-2}});
    mluseDif->AssignWENOStencils(0,"(5,5)",{{-2,-2}});

    mluseDif->UpdateNonLinearWgts(drivPtr->mi, "interior", "two_stage", interior);

    // Finialize the program ====================================================

    // Clear flow vectors
    VecDestroy(&destStokes_sol);
    VecDestroy(&destStokes_g);
    VecDestroy(&destDarcy_sol);
    VecDestroy(&destDarcy_g);

    free(refArrayStokes);
    free(refArrayDarcy);

    DMDAVecRestoreArray(dmu,localu,&lu);
    DMRestoreLocalVector(dmu, &localu); 

    PetscCall(VecDestroy(&globalmesh));
    PetscCall(VecDestroy(&globalu));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    PetscFinalize();

    return 0;
}
