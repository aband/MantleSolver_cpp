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
    cout << "Initialize with vector output ...";

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

	 // Define transport variables
    TransportVariable myH = TransportVariable();
    TransportVariable myC = TransportVariable();

    myC.diffusion = false;
    myH.diffusion = false;

    mycouple->PrepareTransport(myH, myC, InitHD, InitCD);

    // Read in preheating result
    mycouple->ReadVectorTransport(&myH.sol, &myC.sol, "cellH", "cellC", 1);

    mycouple->printCellGrids();

    mycouple->printCellScalar(&myC.sol, "newcellC", 1);
    mycouple->printCellScalar(&myH.sol, "newcellH", 1);

    myH.CreateDefaultReconstruction(mycouple->mi);
    myC.CreateDefaultReconstruction(mycouple->mi);

    myH.Evaluate(mycouple->mi, mycouple->dmu);
    myC.Evaluate(mycouple->mi, mycouple->dmu);
    // ======================================================================
    // Define transport variable	of temperature
    Vec Temp;
	 PetscCall(VecDuplicate(myH.sol, &Temp));
    PetscCall(VecCopy(myH.sol, Temp));

    // Compute porosity with phase package
    mycouple->computePorosity_phase(myH, myC);
    mycouple->adjustEdgePorosity();
    mycouple->assignTempVec(&Temp);

    mycouple->printedgeporosity(1);

    // Diffusion 
    TransportVariable myTempDif = TransportVariable();

    PetscCall(VecDuplicate(Temp, &myTempDif.sol));
    PetscCall(VecCopy(Temp, myTempDif.sol));

    myTempDif.diffusion = true;

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

    myTempDif.CreateReconstruction(mycouple->mi, sizelgx, sizelgy, orderlg, 
			       		                            sizesmx, sizesmy, ordersm,
			    		                               sten_lg_pre, mylinwgts_lg,
			    											    sten_sm_pre, mylinwgts_sm,
			    		                               true, mylinwgts_const); 

    // Advection
    TransportVariable myTempAdv = TransportVariable();
 
    PetscCall(VecDuplicate(Temp, &myTempAdv.sol));
    PetscCall(VecCopy(Temp, myTempAdv.sol));

    myTempAdv.diffusion = false;
    myTempAdv.CreateDefaultReconstruction(mycouple->mi);

    // Evaluate both diffusion and advection
    myTempDif.Evaluate(mycouple->mi, mycouple->dmu);
	 myTempAdv.Evaluate(mycouple->mi, mycouple->dmu);

    Vec latent;
    PetscCall(VecDuplicate(myH.sol, &latent));
	 PetscCall(VecCopy(myH.sol, latent));

    mycouple->setlatentVec(&latent);

    Vec Hsol;
    PetscCall(VecDuplicate(myH.sol, &Hsol));
	 PetscCall(VecCopy(myH.sol, Hsol));

//    VecView(latent, PETSC_VIEWER_STDOUT_WORLD);
/*
    TransportVariable expandporosity = TransportVariable();

	 PetscCall(VecDuplicate(myH.sol, &expandporosity.sol));
    PetscCall(VecCopy(myH.sol, expandporosity.sol));

	 mycouple->AssignPorosityVec(&expandporosity.sol);

    expandporosity.diffusion = false;
//    expandporosity.CreateDefaultReconstruction(mycouple->mi);

    // Reuse former setup

    expandporosity.CreateReconstruction(mycouple->mi, sizelgx, sizelgy, orderlg, 
			       		                            sizesmx, sizesmy, ordersm,
			    		                               sten_lg_pre, mylinwgts_lg,
			    											    sten_sm_pre, mylinwgts_sm,
			    		                               true, mylinwgts_const); 

	 expandporosity.Evaluate(mycouple->mi, mycouple->dmu);
    //mycouple->expandporosity(expandporosity, M+1, N,  M, N, 0);
	 //mycouple->expandporosity(expandporosity, M,   N+1,M, N, (M+1)*N);
*/

    //mycouple->expandporosity();

    // ======================================================================
	 cout << "Initialize Darcy-Stokes solver ... " << endl;
    DarcyStokes ds = DarcyStokes(mycouple->mi, mycouple->myPhase.pp, {0.0});

    // ======================================================================
    int mark = 1;
    int frame = 50;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-frame", &frame, NULL)); 

    //mycouple->computePorosity();

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

    // Only print out three quantities here
    mycouple->printdivmass(mark,"massconsv",ds.StokesVel, ds.DarcyVel);
    mycouple->printphase(mark);
    mycouple->printedgeporosity(mark);

	 /*
    for (int t=0; t<Tmax; t++){

        cout << "Time " << t << endl;
        // Assemble linear system for every time step
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

        if (t%frame ==0){
            mycouple->printdivmass(mark,"massconsv",ds.StokesVel, ds.DarcyVel);
            mycouple->printedgevel(mark, mycouple->effvel, mycouple->phasevel);
				mycouple->printphase(mark);
				mycouple->printedgeporosity(mark);
				mark ++;
        }

        // Set temperature as transport solution and evaluate WENO coefficients
        myC.Evaluate(mycouple->mi,mycouple->dmu);

        // Advect composition
        Vec fluxCAdv; 
		  PetscCall(VecDuplicate(myC.sol, &fluxCAdv));
        myC.advflux_all(mycouple->mi, 1e-5, 
								mycouple->dmu, mycouple->edgegauss, 
								mycouple->effvel, true, 1, &fluxCAdv);       

        // Advect-diffuse enthalpy
        Vec tempfluxDif;
		  PetscCall(VecDuplicate(myH.sol, &tempfluxDif));

        Vec tempfluxAdv;
        PetscCall(VecDuplicate(myH.sol, &tempfluxAdv));
        myTempAdv.advflux_all(mycouple->mi, 1e-5, 
								      mycouple->dmu, mycouple->edgegauss, 
										mycouple->phasevel, true, 0, &tempfluxAdv);

        // Evolution in time
        PetscCall(VecAXPY(myC.sol, -1*dt, fluxCAdv));

        PetscCall(VecAXPY(tempfluxAdv, 8e-8, tempfluxDif)); 
        PetscCall(VecAXPY(Temp, -1*dt, tempfluxAdv));

        // Compute phase variables and adjust porosity
        mycouple->computePorosity_phase(myH, myC);
        //mycouple->adjustEdgePorosity();
        
		  // ===Expand cell centered porosity to edges 
		  cout << "Check " << endl;
        mycouple->AssignPorosityVec(&expandporosity.sol);
        expandporosity.Evaluate(mycouple->mi, mycouple->dmu);
        mycouple->expandporosity(expandporosity, M+1, N,  M, N, 0);
		  mycouple->expandporosity(expandporosity, M,   N+1,M, N, (M+1)*N);
        // ===========================================

        mycouple->assignTempVec(&Temp);

        //mycouple->examineFullPorosity(t);

        PetscCall(VecCopy(Temp, myTempDif.sol));
        PetscCall(VecCopy(Temp, myTempAdv.sol));

        // Evaluate new WENO coefficients
        myC.Evaluate(mycouple->mi, mycouple->dmu);	

        myTempAdv.Evaluate(mycouple->mi, mycouple->dmu);
        myTempDif.Evaluate(mycouple->mi, mycouple->dmu);
    }
*/

    return 1;
}
