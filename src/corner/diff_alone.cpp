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

    double L = 2.0, H = 2.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));

    double addy = 0.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-addy",&addy,NULL));

    //double xstart = -0.5*L, ystart = -1.0001*H - addy;
    double xstart = -1.0, ystart = -1.0;
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
    TransportVariable myC = TransportVariable();

    myC.diffusion = true;
    myH.diffusion = true;

    mycouple->PrepareTransport(myH, myC, InitHD, InitCD);       

    myH.CreateDefaultReconstruction(mycouple->mi);
    //myC.CreateDefaultReconstruction(mycouple->mi);

    int sizelgx = 5;
    int sizelgy = 5;
    int orderlg = 4;

    int sizesmx = 3;
    int sizesmy = 3;
    int ordersm = 2;

    vector<indice> sten_lg_pre = {{-2,-2}};
    vector<indice> sten_sm_pre = {{-2,-2},{-2,0},{0,-2}, {0,0}};

    vector<double> mylinwgts_lg = {0.0};
	 vector<double> mylinwgts_sm = {0.0,0.0,0.0,0.0};
	 double mylinwgts_const = 0.0001;

    myC.CreateReconstruction(mycouple->mi, sizelgx, sizelgy, orderlg, 
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

    mycouple->printCellScalar(&myC.sol, "Init", mark);

    int frame = 1;
//    int count = 0;
    // Euler forwarding
    for (int t=0; t<Tmax; t++){
        myC.Evaluate(mycouple->mi, mycouple->dmu);
        //myC.PrintSample(mycouple->mi);

        Vec fluxC;
        VecDuplicate(myC.sol, &fluxC);
        //myC.cellflux_all(mycouple->mi, 1e-5, mycouple->dmu, mycouple->edgegauss, ds.StokesVel, true, 0, &fluxC);
        myC.difflux_all(mycouple->mi, mycouple->dmu, mycouple->edgegauss, &fluxC);

 //       myC.PrintSample(mycouple->mi, M+1, N, 0, "samplevert.dat", "samplevxvert.dat", "samplevyvert.dat");
 //       myC.PrintSample(mycouple->mi, M, N+1, (M+1)*N, "samplehori.dat", "samplevxhori.dat", "samplevyhori.dat");

        myC.PrintEdgeSample(mycouple->mi, M+1, N, 0, "samplevert.dat", "samplevxvert.dat", "samplevyvert.dat");
        mycouple->printCellScalar(&fluxC, "flux", mark);

        //VecView(fluxH, PETSC_VIEWER_STDOUT_WORLD);
        PetscCall(VecAXPY(myC.sol, -1*dt, fluxC));
		  if (t%frame == 0){
            mycouple->printCellScalar(&myC.sol, "cellC", mark);
				mark ++ ;
		  }

    }
    return 1;
}
