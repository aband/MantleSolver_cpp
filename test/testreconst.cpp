#include <iostream>
#include <petsc.h>
#include "weno.h"
#include "integral.h"
#include "input.h"
#include "flux.h"

#include "func.h"

#include <adolc/adolc.h>

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

void testoutput(const valarray<double> test){
    valarray<double> get = test;
	 /*
     *cout << test.size() << endl;
     *cout << get.size() << endl;
	  */
	 for (auto &i: get){
		  i+=1;
	 }
	 /*
	  *for (auto &i: get){
	  *    cout << i << endl;
	  *}
	  */
}

inline valarray<double> testoutput2(valarray<double> test){
    valarray<double> gettest = test;
    for (auto &i: gettest){
        i = 7+i;
    }
    return gettest;
}

double func(valarray<double>& point, const vector<double>& param){
/*
 *    if (point[0]<0.50){
 *        return sin(point[0])+cos(point[1]);
 *    } else {
 *        return sin(point[0])+cos(point[1])+10;
 *        //return exp(point[0]+point[1]);
 *    }
 *
 */
	 return sin(point[0]+point[1])+cos(exp(point[1]));
    //return point[0]*point[0] + point[1]*point[1];
}

/*
 *double funcX(valarray<double>& target, const vector<double>& param){
 *
 *    return param[0]*param[0]/2;
 *}
 *
 *double funcY(valarray<double>& target, const vector<double>& param){
 *
 *    return param[0]*param[0]/2;
 *}
 *
 *double dfuncX(valarray<double>& target, const vector<double>& param){
 *
 *    return param[0];
 *}
 *
 *double dfuncY(valarray<double>& target, const vector<double>& param){
 *
 *    return param[0];
 *}
 *
 */
int main(int argc, char **argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processors \n",size);CHKERRQ(ierr);

    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // valarray has dynamic memory allocation
/*
 *    int k = 2;
 *    valarray<double> test(k);
 *
 *    test[0] = 1.0;
 *    test[1] = 2.0;
 *
 *    vector< valarray<int> > index { {-1,0,}, {-1,-1}, {0,1}, {1,0} };
 *
 *    wenobasis * work;
 *
 *    work = new wenobasis(test,index,2);
 *
 *    testoutput(test);
 *
 *    valarray<double> duplicate = testoutput2(test);
 *    duplicate = duplicate*duplicate;
 *
 *    // Test for integral function
 *    vector< valarray<double> > corner {{0.0},{1,0.0},{1.5,1.5},{0.0,1}};
 *    valarray<double> center {0.0,0.0};
 *    double h = 1.0;
 *
 *    double result = NumIntegralFace(corner,center,h,func);
 *
 */
//    cout << result << endl;
	 /*
     *cout << duplicate[0] << endl;
     *cout << test[0] << endl;
     *cout << index[0][0] << endl;
	  */
/*
 *    cout << work[0].center[0] << work[0].center[1] << endl;
 *    for (auto & i: index){
 *        cout << i[0] << i[1] << endl;
 *    }
 *
 */

/*
 *    vector<double> Empty;
 *    vector< valarray<double> > corner {{0.0},{1.0/3.0,0.0},{1.0/3.0,1.0/3.0},{0.0,1.0/3.0}};
 *    cout << NumIntegralFace(corner,Empty,{0.0,0.0},1.0,func) << endl;
 *
 */
    // ==========================================================================================================================

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 3, N = 3;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    // Create data management object
    DM    dm;
    Vec   fullmesh;
    const int stencilWidth = 5;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidth, NULL, NULL, &dm);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dm);               CHKERRQ(ierr);
    ierr = DMSetUp(dm);                        CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm, &fullmesh);CHKERRQ(ierr); 

    double L = 3.0, H = 3.0;
    double xstart = -1.0, ystart = -1.0;
    ierr = PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL); CHKERRQ(ierr);

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

    cout << "Created full mesh. To check full mesh, rerun with -printmesh 1 " << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // ==========================================================================================================================

    // Contain defined mesh in vector container.
    // and verify it.
    vector< valarray<double> > mesh;
    
    ReadMeshPortion(dm, &fullmesh, mesh);

    //cout << "Converted c array of local mesh into vector container c++ " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    vector<int> Lorder {3,3};
    vector<int> Sorder {2,2};
    index_set StencilLarge {{-1,-1},{0,-1},{1,-1},
                            {-1, 0},{0, 0},{1, 0},
                            {-1, 1},{0, 1},{1, 1}};

    vector<index_set> StencilSmall;
    StencilSmall.push_back({{-1,-1},{ 0,-1},{ 0, 0},{-1, 0}});
    StencilSmall.push_back({{ 0,-1},{ 1,-1},{ 1, 0},{ 0, 0}});
    StencilSmall.push_back({{ 0, 0},{ 1, 0},{ 1, 1},{ 0, 1}});
    StencilSmall.push_back({{-1, 0},{ 0, 0},{ 0, 1},{-1, 1 }});

/*
 *    typedef WenoReconst*  wrPtr;
 *
 *    wrPtr * wr_23 = new wrPtr[StencilNum];
 *
 *    for (int s=0; s<StencilNum; s++){
 *        int j=s/(M+2);
 *        int i=s%(M+2);
 *        point_index target {i-1+wm->ghost, j-1+wm->ghost};
 *        wr_23[s] = new WenoReconst(target, wm, StencilLarge, Lorder, StencilSmall, Sorder);
 *        wr_23[s]->CreateCoefficients();
 *    }
 *
 *    // Test with new weno reconst object
 *    solution u_reconst = wr_23[(N/2+1)*(M+2)+(M/2+1)]->PointReconstruction(wm, {0.5,0.5});
 *
 *    point ref {0.5,0.5};
 *
 *    cout << endl << "The exact value of the given point (0.5,0.5) is : " << func(ref,{0.0}) << endl;
 *
 *    cout << endl << "The reconstruction of the point (0.5,0.5) is : " << u_reconst << endl;
 *
 */

    // ==========================================================================================================================

    DM dmu;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, 2, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu);               CHKERRQ(ierr);
    ierr = DMSetUp(dmu);                        CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    //SimpleInitialValue(dm,dmu,&fullmesh,&globalu,func);
    // Initialize with oblique data for Burgers equation 
    ObliqueBurgers(dm,dmu,&fullmesh,&globalu,func);

    Vec localu; 
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu);

    // It can be changed later to not be double
    solution  ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    /*
     *Initialize code by calling WenoMesh
     */
    const WenoMesh * wm = new WenoMesh(M,N,stencilWidth,mesh,lu);

    PetscInt       xs,ys,xm,ym;
    ierr = DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

    // Create array of weno prepare object
    // First determine how many stencils we need locally
    int StencilNum = (M+2)*(N+2);

	 // Benchmark calculation
    WenoReconst * wr_23;

    point_index target_index = {N/2+wm->ghost ,M/2+wm->ghost};

    wr_23 = new WenoReconst(target_index,wm,StencilLarge,Lorder,StencilSmall,Sorder);
    wr_23->CreateNewWeights();
    //wr_23->CreateCoefficients();
    wr_23->CheckWeights();
    
	 point target_point = {0.5, 0.5}; 

	 printf("\nThe target point is (%.2f, %.2f), the exact value is %.12f \n \n", target_point[0], target_point[1], func(target_point,{0.0}));

	 // Define an instance for 3,2 reconstruction

    solution u = wr_23->PointReconstruction(wm, {0.5,0.5});

	 printf("The reconstructed value at the given point is %.12f \n", u);

	 printf("\nThe error measured at this point : %.12f \n\n",abs(u-func(target_point,{0.0})));

    delete wr_23;

    DMDAVecRestoreArray(dmu,localu,&lu);
 
    // ==========================================================================

    // Destroy Vectors
    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    DMDestroy(&dm);

    return 0;
}
