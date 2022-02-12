#include <iostream>
#include <petsc.h>
#include "../include/weno.h"
#include "../include/integral.h"
#include "../include/input.h"

#include <adolc/adolc.h>

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

adouble power(adouble x, int n)
{
    adouble z = 1;
    if (n>0) {
        int nh = n/2;
        z = power(x,nh);
        z *= z;
        if (2*nh !=n) {z *= x;}
        return z;
    } else {
        if (n==0) {
            return z;
        } else {
            return 1.0/power(x,-n);
        }
    }
}

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

double func(valarray<double>& point){
    //return point[0]*point[0];
    return sin(point[0])+cos(point[1]);
}

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

    double L = 1.0, H = 1.0;
    double xstart = 0.0, ystart = 0.0;
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

/*
 *    for (auto & i: mesh){
 *        cout << "(" << i[0] << "," << i[1] << ")" << endl;
 *    }
 *
 */

//    WenoMesh wm = WenoMesh(M,N,stencilWidth,mesh);

    cout << "Convert c array of local mesh into vector container c++ " << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // ==========================================================================================================================

    
    // WENO2 stencil and corresponding reconstruction.
    vector< valarray<int> > stencilindex {{-1,-1},{0,-1},{0,0},{-1,0}};
    vector<int> order {2,2};
    vector<int> range {stencilWidth-1,stencilWidth+M+2,
                       stencilWidth-1,stencilWidth+N+2};

    WenoBasisCoeff * wbc_weno2;
    wbc_weno2 = new WenoBasisCoeff(M,N,stencilWidth,mesh);

    wbc_weno2->SetStencilInfo(stencilindex,order,range);

    wbc_weno2->GetStencil();

    wbc_weno2->CreateWenoBasisCoeff();

    // Destroy Vectors
    VecDestroy(&fullmesh);
    DMDestroy(&dm);

    return 0;
}
