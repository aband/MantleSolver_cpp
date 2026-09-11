#include "couple.h"

int main(int argc, char **argv){

    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Input mesh parameter =========================================================
    int M=4, N=4;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL));

    double L = 0.5, H = 0.5;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));

    double addy = 0.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-addy",&addy,NULL));

    //double xstart = -0.5*L, ystart = -1.0001*H - addy;
    double xstart = 0.0, ystart = -1.000*H;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    int stencilWidthMesh = 5;
    int stencilWidthU = 3;

    int physicsScale = 0;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-scale", &physicsScale, NULL));

    int meshType = 0; 
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshType,NULL));

    int maxIter = 100; 
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL));    

    double tolUzawa = 10e-14; 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL)); 

    double Tmax = 20; // Stop at the first step 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tmax", &Tmax, NULL)); 

    double dt = 0;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));

    int showPhase = 0;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-showphase", &showPhase, NULL)); 

    int withUnit = 0;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-unit", &withUnit, NULL));

    int interval = 1;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-interval", &interval, NULL));

    // Read vector and assign as initial value
    cout << "Initialize with vector output ...";

    couple * mycouple = new couple(); 

    // Initialize coupling variables
    mycouple->withUnit = withUnit;
    mycouple->CreatePhase();
    mycouple->ShowPhase();
    mycouple->CreateMesh(M, N, L, H, xstart, ystart, 
                         stencilWidthMesh, stencilWidthU,
                         physicsScale, meshType);

	 // Define transport variables
    TransportVariable myH = TransportVariable();
    TransportVariable myC = TransportVariable();

    myC.diffusion = false;
    myH.diffusion = false;

    mycouple->PrepareTransport(myH, myC, InitHD, InitCD);

    mycouple->ReadVectorTransport(&myH.sol, &myC.sol, "cellH", "cellC", 1);

    mycouple->printCellScalar(&myC.sol, "newcellC", 1);
    mycouple->printCellScalar(&myH.sol, "newcellH", 1);

    myH.CreateDefaultReconstruction(mycouple->mi);
    myC.CreateDefaultReconstruction(mycouple->mi);

    myH.Evaluate(mycouple->mi, mycouple->dmu);
    myC.Evaluate(mycouple->mi, mycouple->dmu);

    // Initialize porosity according to phase calculation
    mycouple->computePorosity_phase(myH, myC);
    mycouple->printedgeporosity(1);
    mycouple->printphase(1);

    // Compute velocity
	 cout << "Initialize Darcy-Stokes solver ... " << endl;
    DarcyStokes ds = DarcyStokes(mycouple->mi, mycouple->myPhase.pp, {0.0});

    ds.Assemble(mycouple->mi, 
                mycouple->edgeporo,
                mycouple->cellporo,
                mycouple->average_poro,
                0.0,
                mycouple->myPhase.pp);

    ds.CreateCoupledSystem();

    ds.Solve(maxIter, tolUzawa);

    ds.ReconstructEdgeVel(mycouple->edgegauss, mycouple->mi);

    mycouple->calculatePhaseVel(ds.StokesVel, ds.DarcyVel);

//    mycouple->printedgevel(1, mycouple->effvel, mycouple->phasevel);
    mycouple->printedgevel(1, ds.StokesVel, ds.DarcyVel);

    Vec sp, dp;
    ds.PreparePressure(&sp, &dp);

    // Define additional diffusion and latent heat transport term
    TransportVariable myHdiff = TransportVariable();
    TransportVariable myHLatent = TransportVariable();

    myHLatent.diffusion = false;

	 myHdiff.diffusion = true;

    // Copy solution
    PetscCall(VecDuplicate(myH.sol, &myHdiff.sol));
    PetscCall(VecCopy(myH.sol, myHdiff.sol));

    PetscCall(VecDuplicate(myH.sol, &myHLatent.sol));
    PetscCall(VecCopy(myH.sol, myHLatent.sol));

    myHLatent.CreateDefaultReconstruction(mycouple->mi);

    // Diffusion weno stencils
    int sizelgx = 3;
    int sizelgy = 3;
    int orderlg = 2;

    int sizesmx = 2;
    int sizesmy = 2;
    int ordersm = 1;

    vector<indice> sten_lg_pre = {{-1,-1}};
    vector<indice> sten_sm_pre = {{0,0}};

    vector<double> mylinwgts_lg = {0.0};
	 vector<double> mylinwgts_sm = {0.0,};
	 double mylinwgts_const = 1;

    myHdiff.CreateReconstruction(mycouple->mi, sizelgx, sizelgy, orderlg, 
			     		                             sizesmx, sizesmy, ordersm,
					                                sten_lg_pre, mylinwgts_lg,
														     sten_sm_pre, mylinwgts_sm,
					                                true, mylinwgts_const); 

    int mark = 2;

    Vec fluxC;
    PetscCall(VecDuplicate(myC.sol, &fluxC));

    Vec fluxH;
    PetscCall(VecDuplicate(myH.sol, &fluxH));

    Vec temp;
    PetscCall(VecDuplicate(myH.sol, &temp));

    // Coupled time stepping
    for (int t=0; t<Tmax; t++){
        myC.advflux_all(mycouple->mi, 1e-5, mycouple->dmu, mycouple->edgegauss, 
								mycouple->effvel, true, 0, &fluxC);

        myH.advflux_all(mycouple->mi, 1e-5, mycouple->dmu, mycouple->edgegauss, 
								mycouple->phasevel, true, 0, &fluxH);

//        myH.flux

        mycouple->assignTempVec(&temp); 

        myC.Evaluate(mycouple->mi, mycouple->dmu);
        myH.Evaluate(mycouple->mi, mycouple->dmu);

	 }

    return 1;
}
