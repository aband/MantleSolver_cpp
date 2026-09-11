// Standardized driver function that activates all the function

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

    // ==============================================================================
    couple * mycouple = new couple(); 

    // Initialize coupling variables
    mycouple->withUnit = withUnit;
    mycouple->CreatePhase();
    mycouple->ShowPhase();
    mycouple->CreateMesh(M, N, L, H, xstart, ystart, 
                         stencilWidthMesh, stencilWidthU,
                         physicsScale, meshType);

    mycouple->printGaussPoints();
    mycouple->printCellGrids();

    mycouple->computePorosity();

    //mycouple->printedgeporosity(1);

    // Initialize darcy stokes solver
    DarcyStokes ds = DarcyStokes(mycouple->mi, mycouple->myPhase.pp, {0.0});

    ds.Assemble(mycouple->mi, 
                mycouple->edgeporo,
                mycouple->cellporo,
                mycouple->average_poro,
                0.0,
                mycouple->myPhase.pp);

    //ds.showMatrix();
    //ds.printBndryAll();

    ds.CreateCoupledSystem();

    ds.Solve(maxIter, tolUzawa);

    ds.ReconstructEdgeVel(mycouple->edgegauss, mycouple->mi);

    mycouple->printedgevel(1, ds.StokesVel, ds.DarcyVel);

    cout << "Transport is enabled here." << endl;

    TransportVariable myH = TransportVariable();
    TransportVariable myCAdv = TransportVariable();

    myCAdv.diffusion = false;
    myH.diffusion = false;

    mycouple->PrepareTransport(myH, myCAdv, InitHD, InitCD);       

    // Advection weno stencils
    myH.CreateDefaultReconstruction(mycouple->mi);
    myCAdv.CreateDefaultReconstruction(mycouple->mi);

    // =====================================================

    TransportVariable myCDif = TransportVariable();

    PetscCall(VecDuplicate(myCAdv.sol, &myCDif.sol));
    PetscCall(VecCopy(myCAdv.sol, myCDif.sol));

    myCDif.diffusion = true;
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

    myCDif.CreateReconstruction(mycouple->mi, sizelgx, sizelgy, orderlg, 
			     		                            sizesmx, sizesmy, ordersm,
					                               sten_lg_pre, mylinwgts_lg,
														    sten_sm_pre, mylinwgts_sm,
					                               true, mylinwgts_const); 

//    myH.Evaluate(mycouple->mi,mycouple->dmu);        
//    myC.Evaluate(mycouple->mi,mycouple->dmu);

//    myH.Print(mycouple->mi, GetFilename("H", 1));
//    myC.Print(mycouple->mi, GetFilename("C", 1));

//    Vec fluxH;
//    VecDuplicate(myH.sol,&fluxH);
//    myH.cellflux_all(mycouple->mi, 1e-4, mycouple->dmu, mycouple->edgegauss, ds.StokesVel, false, 0, &fluxH);

//    mycouple->printCellScalar(&myH.sol, "cellH", 1);
//    mycouple->printCellScalar(&myH.sol, "cellH", 2);

    int mark = 1;

    int frame = 50;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-frame", &frame, NULL)); 

    // Euler forwarding
    for (int t=0; t<Tmax; t++){
        myCAdv.Evaluate(mycouple->mi, mycouple->dmu);
        myCDif.Evaluate(mycouple->mi, mycouple->dmu);
        //myC.PrintSample(mycouple->mi);

        Vec fluxCAdv;
        PetscCall(VecDuplicate(myCAdv.sol, &fluxCAdv));
        myCAdv.advflux_all(mycouple->mi, 1e-5, mycouple->dmu, mycouple->edgegauss, 
								ds.StokesVel, true, 0, &fluxCAdv);

        Vec fluxCDiv;
        PetscCall(VecDuplicate(myCDif.sol, &fluxCDiv));
        myCDif.difflux_all(mycouple->mi, mycouple->dmu, mycouple->edgegauss, &fluxCDiv);

        //PetscCall(VecAXPY(fluxCAdv, 5e-9, fluxCDiv));
        PetscCall(VecAXPY(fluxCAdv, 8e-8, fluxCDiv));

        //VecView(fluxH, PETSC_VIEWER_STDOUT_WORLD);
        PetscCall(VecAXPY(myCAdv.sol, -1*dt, fluxCAdv));

//		  if (t%frame == 0){
//            mycouple->printCellScalar(&myCAdv.sol, "cellH", mark);
//				mark ++ ;
//		  }

        PetscCall(VecCopy(myCAdv.sol, myCDif.sol));

    }

    cout << "Pre heating finished ... " << endl;
    //mycouple->printCellScalar(&myCAdv.sol, "cellH", 1);
    //mycouple->printCellScalar(&myH.sol, "cellC", 1);

    myCAdv.Evaluate(mycouple->mi, mycouple->dmu);
    myH.Evaluate(mycouple->mi, mycouple->dmu);

    mycouple->computePorosity_phase(myCAdv, myH);
    mycouple->printedgeporosity(1);
    mycouple->printphase(1);

    // Assemble darcy-stokes system
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

    mycouple->printedgevel(1, ds.StokesVel, ds.DarcyVel);

    mycouple->printCellScalar(&myH.sol, "cellH", 1);
    mycouple->printCellScalar(&myH.sol, "cellH", 2);

    return 1;
}
