// Coupled porosity evolving file
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
	 cout <<"Initialize coupling variables ..." << endl; 
    mycouple->withUnit = withUnit;
    mycouple->CreatePhase();
    mycouple->ShowPhase();
    mycouple->CreateMesh(M, N, L, H, xstart, ystart, 
                         stencilWidthMesh, stencilWidthU,
                         physicsScale, meshType);

    mycouple->printGaussPoints();
    mycouple->printCellGrids();

    cout << "Define transport variables ... " << endl;
    TransportVariable myC = TransportVariable();
    TransportVariable myHAdv = TransportVariable();

    myHAdv.diffusion = false;
    myC.diffusion = false;

    mycouple->PrepareTransport(myHAdv, myC, InitHD, InitCD);       

    // Advection weno stencils
    myC.CreateDefaultReconstruction(mycouple->mi);
    myHAdv.CreateDefaultReconstruction(mycouple->mi);

    // =====================================================

    TransportVariable myHDif = TransportVariable();

    PetscCall(VecDuplicate(myHAdv.sol, &myHDif.sol));
    PetscCall(VecCopy(myHAdv.sol, myHDif.sol));

    myHDif.diffusion = true;
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

    myHDif.CreateReconstruction(mycouple->mi, sizelgx, sizelgy, orderlg, 
			     		                            sizesmx, sizesmy, ordersm,
					                               sten_lg_pre, mylinwgts_lg,
														    sten_sm_pre, mylinwgts_sm,
					                               true, mylinwgts_const); 

    myC.Evaluate(mycouple->mi, mycouple->dmu);
    myHAdv.Evaluate(mycouple->mi, mycouple->dmu);
    myHDif.Evaluate(mycouple->mi, mycouple->dmu);

//    cout << "Compute initial porosity distribution ... " << endl;
//    mycouple->computePorosity_phase(myH, myCAdv);

    cout << "Initialize Darcy-Stokes solver ... " << endl;
    DarcyStokes ds = DarcyStokes(mycouple->mi, mycouple->myPhase.pp, {0.0});

    cout << "Start time evolution ... " << endl;

    int mark = 1;
    int frame = 50;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-frame",&frame,NULL));

	 // Solve static Stokes problem
    mycouple->computePorosity();

    ds.Assemble(mycouple->mi, 
                mycouple->edgeporo,
                mycouple->cellporo,
                mycouple->average_poro,
                0.0,
                mycouple->myPhase.pp);

    ds.CreateCoupledSystem();

    ds.Solve(maxIter, tolUzawa);

    ds.ReconstructEdgeVel(mycouple->edgegauss, mycouple->mi);

    mycouple->printedgevel(1, ds.StokesVel, ds.DarcyVel);

// ===========================================================================
    int mark = 1;

    int frame = 50;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-frame",&frame,NULL));

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


//        PetscCall(VecView(fluxCAdv, PETSC_VIEWER_STDOUT_WORLD));
//        PetscCall(VecView(fluxCDiv, PETSC_VIEWER_STDOUT_WORLD));

        PetscCall(VecAXPY(fluxCAdv, 0.8e-7, fluxCDiv));

        //VecView(fluxH, PETSC_VIEWER_STDOUT_WORLD);
        PetscCall(VecAXPY(myCAdv.sol, -1*dt, fluxCAdv));

		  if (t%frame == 0){
            mycouple->printCellScalar(&myCAdv.sol, "cellC", mark);
				mark ++ ;
		  }

        PetscCall(VecCopy(myCAdv.sol, myCDif.sol));

    }


	 /*
    // Precompute the system
    for (int t1=0; t1<100; t1++){
        myHAdv.Evaluate(mycouple->mi, mycouple->dmu);
        myHDif.Evaluate(mycouple->mi, mycouple->dmu);
 
        Vec fluxHAdv;
        PetscCall(VecDuplicate(myHAdv.sol, &fluxHAdv));
        myHAdv.advflux_all(mycouple->mi, 1e-5, mycouple->dmu, mycouple->edgegauss, 
								ds.StokesVel, true, 0, &fluxHAdv);

        Vec fluxHDiv;
        PetscCall(VecDuplicate(myHDif.sol, &fluxHDiv));
        myHDif.difflux_all(mycouple->mi, mycouple->dmu, mycouple->edgegauss, &fluxHDiv);

        PetscCall(VecAXPY(fluxHAdv, 0.8e-7, fluxHDiv));

        PetscCall(VecAXPY(myHAdv.sol, -1*100, fluxHAdv));

        PetscCall(VecCopy(myHAdv.sol, myHDif.sol));

//        mycouple->printedgeporosity(mark);
	     mycouple->printCellScalar(&myHAdv.sol, "cellH", mark);
		  mark ++ ;
	
	 }

	 for (int t=0; t<Tmax; t++){

        // Evaluate porosity distribution
        mycouple->computePorosity_phase(myHAdv, myC);

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

//        mycouple->printedgevel(mark, mycouple->effvel, mycouple->phasevel);


        // Initialize flux
        Vec fluxHAdv;
        PetscCall(VecDuplicate(myHAdv.sol, &fluxHAdv));
        myHAdv.advflux_all(mycouple->mi, 1e-5, mycouple->dmu, mycouple->edgegauss, 
								mycouple->phasevel, true, 0, &fluxHAdv);

        Vec fluxHDiv;
        PetscCall(VecDuplicate(myHDif.sol, &fluxHDiv));
        myHDif.difflux_all(mycouple->mi, mycouple->dmu, mycouple->edgegauss, &fluxHDiv);

        Vec fluxC;
        PetscCall(VecDuplicate(myC.sol, &fluxC));
        myHAdv.advflux_all(mycouple->mi, 1e-5, mycouple->dmu, mycouple->edgegauss, 
								mycouple->effvel, true, 1, &fluxC);

        PetscCall(VecAXPY(fluxHAdv, 0.8e-7, fluxHDiv));

        // Time stepping
        PetscCall(VecAXPY(myHAdv.sol, -1*dt, fluxHAdv));

        PetscCall(VecAXPY(myC.sol, -1*dt, fluxC));

//		  if (t%frame == 0){
//            mycouple->printCellScalar(&myHAdv.sol, "cellH", mark);
//				mark ++ ;
//		  }

        PetscCall(VecCopy(myHAdv.sol, myHDif.sol));

	 }	
*/
    return 1;
}
