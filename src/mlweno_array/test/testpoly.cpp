#include "petsc.h"
#include "input.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

#include "tensorstencilpoly.h"
#include "reconstruction.h"
#include <chrono>

#include "error.h"

double func(const vertex& point,
            const vector<double>& param){

/*
    if (point[0] < 0.5) {

    return sin(point[0])*cos(point[1]);

    } else {

    return sin(point[0])*cos(point[1]) +0.5;

    }
*/

    double HD = 2.9+2.5*0.2-2.5*pow(point[1]+0.2,1);

    if (point[1] < -0.20){HD = 2.9 + 2.5*0.20;}

    return HD;
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

    L = 0.1, H = 0.4;
    xstart = -0.5*L, ystart = -1.0001*H;
 
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

/*
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
*/

    // Set full stencils
/*
    cout << M << "  " << N << endl;
    vector<tensorstencilpoly> stenlg;
    int Nlg = N-4;
    int Mlg = M-4;
    stenlg.resize(Mlg*Nlg);
	 int lgx = 5;
	 int lgy = 5;
	 int lgr = 4;

    vector<tensorstencilpoly> stensm;
	 int Nsm = N-2;
	 int Msm = M-2;
    stensm.resize(Msm*Nsm);
	 int smx = 3;
    int smy = 3;
	 int smr = 2;

    vector<indice> sten_lg_pre = {{-2,-2}};
    vector<indice> sten_sm_pre = {{-2,-2},{-2, 0},{0 ,-2},{0,0}};
*/

    cout << M << "  " << N << endl;
    vector<tensorstencilpoly> stenlg;
    int Nlg = N-2;
    int Mlg = M-2;
    stenlg.resize(Mlg*Nlg);
	 int lgx = 3;
	 int lgy = 3;
	 int lgr = 2;

    vector<tensorstencilpoly> stensm;
	 int Nsm = N-1;
	 int Msm = M-1;
    stensm.resize(Msm*Nsm);
	 int smx = 2;
    int smy = 2;
	 int smr = 1;

    vector<indice> sten_lg_pre = {{-1,-1}};
    vector<indice> sten_sm_pre = {{-1,-1},{-1, 0},{0 ,-1},{0,0}};

	 auto start = std::chrono::steady_clock::now();
    for (int j=0; j<Nlg; j++){
    for (int i=0; i<Mlg; i++){
        int s = j*(Mlg)+i;
        stenlg.at(s) = tensorstencilpoly(lgr);
        stenlg.at(s).setCoef(mi,i,j);
		  stenlg.at(s).setSigma();
		  stenlg.at(s).startx = i;
		  stenlg.at(s).starty = j;
    }}

    for (int j=0; j<Nsm; j++){
    for (int i=0; i<Msm; i++){
        int s = j*(Msm)+i;
        stensm.at(s) = tensorstencilpoly(smr);
        stensm.at(s).setCoef(mi,i,j);
		  stensm.at(s).setSigma();
		  stensm.at(s).startx = i;
		  stensm.at(s).starty = j;
    }}

	 auto end = std::chrono::steady_clock::now();
	 auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    cout << "Time Test: " << duration.count() << " ms." << endl;

    stensm.at(0).printCoef();

    // Set all the reconstruction
    vector<reconstruction> my_recon;
    my_recon.resize(M*N);

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
        my_recon.at(s) = reconstruction();

        my_recon.at(s).use_sten_const = 0;

        my_recon.at(s).init(smx,smy,lgx,lgy,smr,lgr,sten_lg_pre, sten_sm_pre, mi,{i,j});
        
    }}

    // Reconstruction test
    Vec globalvec, localvec;
    double ** locvals;

    PetscCall(DMCreateGlobalVector(dmu, &globalvec));

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalvec, {0.0}, func);
    // VecView(globalvec, PETSC_VIEWER_STDOUT_WORLD);
    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvec)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalvec, INSERT_VALUES, localvec));
    PetscCall(DMGlobalToLocalEnd(dmu, globalvec, INSERT_VALUES, localvec));

    PetscCall(DMDAVecGetArray(dmu, localvec, &locvals));

    // Compute sigma 
    vector<double> sigma_lg;
    sigma_lg.resize(stenlg.size());

    vector<double> sigma_sm;
    sigma_sm.resize(stensm.size());

    for (int j=0; j<Nlg; j++){
    for (int i=0; i<Mlg; i++){
        int s = j*(Mlg)+i;
        sigma_lg.at(s) = stenlg.at(s).sigma(locvals); 
    }}   

    for (int j=0; j<Nsm; j++){
    for (int i=0; i<Msm; i++){
        int s = j*(Msm)+i;
        sigma_sm.at(s) = stensm.at(s).sigma(locvals);
    }}

/*
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
        cout << my_recon.at(s).sten_lg.size() + my_recon.at(s).sten_sm.size() + my_recon.at(s).use_sten_const << endl;
    }}
*/
   
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
        my_recon.at(s).extractsigma(sigma_lg, sigma_sm);
        my_recon.at(s).setWgts(1.0/(double)M/(double)N);
		  //my_recon.at(s).printinfo();
        //my_recon.at(s).efforder();
        cout << my_recon.at(s).efforder() << "  ";
    }cout << endl;}

    int midM = M/2;
    int midN = N/2;

cout<< "At the middle cell : " << midM << " , " << midN << endl;
    my_recon.at(midN*M+midM).printinfo();
    my_recon.at(midN*M+midM).printmoreinfo(locvals, stenlg, stensm);


    printexactsol(mi, 0, func, 1, true, {0.0});

    vector<vertex> sample = {{-1+1e-3,-1+1e-3},
                             { 0     ,-1+1e-3},
                             { 1-1e-3,-1+1e-3},
									  {-1+1e-3, 0},
                             { 0     , 0},
                             { 1-1e-3, 0},
									  {-1+1e-3, 1-1e-3},
                             { 0     , 1-1e-3},
                             { 1-1e-3, 1-1e-3}};

    vector<vertex> mapped; 
    mapped.resize(sample.size());

    // Evalutation at sample points
    for (int j=0; j<N; j++){

        for (int i=0; i<M; i++){

            // Extract corners of the selected cell
            vertexSet corners = extractCorners(mi, {i,j});

            for (int g=0; g<sample.size(); g++){
                mapped.at(g) = GaussMapPointsFace(sample.at(g), corners);
            }

            my_recon.at(j*M+i).eval(locvals, mapped, stenlg, stensm);
            //if (j==midN && i==midM){
            //    my_recon.at(j*M+i).eval(locvals, mapped, sten5, sten3);
            //}
        }
    }

    printreconsol(my_recon, M ,N, 1);

    // =================================================================
    DMDAVecRestoreArray(dmu,localvec,&locvals);
    DMRestoreLocalVector(dmu, &localvec); 

    PetscCall(VecDestroy(&globalvec));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    return 1;
}
