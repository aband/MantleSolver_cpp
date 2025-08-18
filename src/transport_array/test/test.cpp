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

#include "advectiveflux.h"

#include "rk.h"

double init(const vertex& point,
            const vector<double>& param){
/*
    if (point[0] < 0.5) {

    return sin(point[0])*cos(point[1]);

    } else {

    return sin(point[0])*cos(point[1]) + 0.5;

    }
*/

return 0.0;

}

int main(int argc, char ** argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    int M = 20, N = 20;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    double dt = 0.2*1.0/(double)M;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-dt", &dt, NULL));

    int Nt = 10;
    ierr = PetscOptionsGetInt(NULL,NULL,"-Nt",&Nt,NULL);CHKERRQ(ierr);

    //Nt *= M;

    double CFL = dt/(1.0/(double)M);

    cout << "dt, dh = " << dt << " , " << 1.0/(double)M << ". " << "CFL number is : " << CFL << ". The final time is : " << 
    Nt * dt << ". " << endl;

    double L = 1, H = 1;
    //double xstart = -L/2, ystart = -H/2;
    double xstart = 0.0, ystart = 0.0;
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

    // ====================================================
    Vec globalvec;

    PetscCall(DMCreateGlobalVector(dmu, &globalvec));

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalvec, {0.0}, init);

    // ====================================================
    cout <<"The total cells : " << M << ", " << N << endl;
    vector<tensorstencilpoly> sten5;
   sten5.resize((M-4)*(N-4));

    vector<tensorstencilpoly> sten3;
    sten3.resize((M-2)*(N-2));

    vector<tensorstencilpoly> sten2;
    sten2.resize((M-1)*(N-1));

//	 auto start = std::chrono::steady_clock::now();
//    for (int j=0; j<N-4; j++){
//    for (int i=0; i<M-4; i++){
//        int s = j*(M-4)+i;
//        sten5.at(s) = tensorstencilpoly(4);
//        sten5.at(s).setCoef(mi,i,j);
//		  sten5.at(s).setSigma();
//		  sten5.at(s).startx = i;
//		  sten5.at(s).starty = j;
//    }}

    for (int j=0; j<N-2; j++){
    for (int i=0; i<M-2; i++){
        int s = j*(M-2)+i;
        sten3.at(s) = tensorstencilpoly(2);
        sten3.at(s).setCoef(mi,i,j);
		  sten3.at(s).setSigma();
		  sten3.at(s).startx = i;
		  sten3.at(s).starty = j;
    }}

    for (int j=0; j<N-1; j++){
    for (int i=0; i<M-1; i++){
        int s = j*(M-1) + i;
        sten2.at(s) = tensorstencilpoly(1);
        sten2.at(s).setCoef(mi,i,j);
		  sten2.at(s).setSigma();
		  sten2.at(s).startx = i;
		  sten2.at(s).starty = j;
    }}

//	 auto end = std::chrono::steady_clock::now();
//	 auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
//    cout << "Time Test: " << duration.count() << " ms." << endl;

    // Prepare for stencils (5,3) reconstruction
    //vector<indice> sten_lg_pre = {{-2,-2}};
    //vector<indice> sten_sm_pre = {{-2,-2},{-2, 0},{0 ,-2},{0,0}, {-1,-1}, {0,-1}};

    // Prepare for stencils (3,2) reconstruction
    vector<indice> sten_lg_pre = {{-1,-1}};
    vector<indice> sten_sm_pre = {{-1,-1}, {-1,0}, {0,-1}, {0,0}};

    vector<reconstruction> my_recon;
    my_recon.resize(M*N);

    // Initializing reconstrucitons for each cell
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
 
        my_recon.at(s) = reconstruction();
        if (i==0){
            my_recon.at(s).use_sten_const = 1;
		  }else {
            my_recon.at(s).use_sten_const = 0;
        }
        my_recon.at(s).init(2,2,3,3,1,2,sten_lg_pre, sten_sm_pre, mi,{i,j});
        //my_recon.at(s).init(3,3,5,5,2,4,sten_lg_pre, sten_sm_pre, mi, {i,j});
    }}

    printexactsol(mi, 0, init, 1, true, {0.0});

//    rk1(dt, Nt, &globalvec, mi, dmu, dmMesh, my_recon, sten3, sten2);
    rk2(dt, Nt, &globalvec, mi, dmu, dmMesh, my_recon, sten3, sten2);

    Vec localvec; 
    double ** locvals;

    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvec)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalvec, INSERT_VALUES, localvec));
    PetscCall(DMGlobalToLocalEnd(dmu, globalvec, INSERT_VALUES, localvec));

    PetscCall(DMDAVecGetArray(dmu, localvec, &locvals));

//    printreconsol(my_recon, M, N, 1, sten3, sten2, mi, locvals);
    printreconsol2(my_recon, M, N, 1, sten3, sten2, mi, locvals);
//    printreconsol2(my_recon, M, N, 1, sten5, sten3, mi, locvals);


    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
				int s = j*M+i;
        cout << my_recon.at(s).efforder() << "  ";
    }cout << endl;}

    DMDAVecRestoreArray(dmu,localvec,&locvals);
    DMRestoreLocalVector(dmu, &localvec); 

    // =======================================
    PetscCall(VecDestroy(&globalvec));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    return 1;
}
