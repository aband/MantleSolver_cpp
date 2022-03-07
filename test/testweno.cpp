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

double func(valarray<double>& point, const vector<double>& param){
	 if (point[0]<0.50-param[0]){
		  return sin(point[0])+cos(point[1]);
	 } else {
		  return sin(point[0])+cos(point[1])+10;
		  //return exp(point[0]+point[1]);
	 }

    //return sin(point[0])+cos(point[1]);
}

double funcX(valarray<double>& target, const vector<double>& param){

    return param[0]*param[0]/2;
}

double funcY(valarray<double>& target, const vector<double>& param){

    return param[0]*param[0]/2;
}

double dfuncX(valarray<double>& target, const vector<double>& param){

    return param[0];
}

double dfuncY(valarray<double>& target, const vector<double>& param){

    return param[0];
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

    //cout << "Converted c array of local mesh into vector container c++ " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // ==========================================================================================================================

    DM dmu;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, 2, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu);               CHKERRQ(ierr);
    ierr = DMSetUp(dmu);                        CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    //SimpleInitialValue(dm,dmu,&fullmesh,&globalu,func);
    // Initialize with oblique data for Burgers equation
    SimpleInitialValue(dm,dmu,&fullmesh,&globalu,Initial_Condition);

    Vec localu; 
    DMCreateLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu);

    // It can be changed later to not be double
    solution  ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    /*
     *Initialize code by calling WenoMesh
     */
    const WenoMesh * wm = new WenoMesh(M,N,stencilWidth,mesh,lu);

    /*
     *Test the center point of the entire domain because it is fixed.
     *In order to get a fixed point, M and N should be odd numbers.
     */
/*
 *    point_index input_target_cell {stencilWidth+M/2,stencilWidth+N/2};
 *    vector<int> Lorder {3,3};
 *    vector<int> Sorder {2,2};
 *    index_set StencilLarge {{-1,-1},{0,-1},{1,-1},
 *                            {-1, 0},{0, 0},{1, 0},
 *                            {-1, 1},{0, 1},{1, 1}};
 *
 *    vector<index_set> StencilSmall;
 *    StencilSmall.push_back({{-1,-1},{ 0,-1},{ 0, 0},{-1, 0}});
 *    StencilSmall.push_back({{ 0,-1},{ 1,-1},{ 1, 0},{ 0, 0}});
 *    StencilSmall.push_back({{ 0, 0},{ 1, 0},{ 1, 1},{ 0, 1}});
 *    StencilSmall.push_back({{-1, 0},{ 0, 0},{ 0, 1},{-1, 1 }});
 *
 *    WenoPrepare * wp = new WenoPrepare(StencilLarge, input_target_cell, Lorder);
 *    wp->SetUpStencil(wm);
 *    wp->CreateBasisCoeff(wm);
 *
 */
/*
 *    for (auto & s: StencilSmall){
 *        WenoPrepare * swp = new WenoPrepare(s,input_target_cell,Sorder);
 *        swp->SetUpStencil(wm);
 *        swp->CreateBasisCoeff(wm);
 *        swp->CreateSmoothnessIndicator(wm,2.0,2.0);
 *        cout << swp->sigma << endl;
 *    }
 *
 */
/*
 *    printf("\nThe target point is (%.2f, %.2f), the exact value is %.12f \n \n", wp->center[0], wp->center[1], func(wp->center,{0.0}));
 *
 *    // Define an instance for 3,2 reconstruction
 *    solution u = WenoPointReconst(StencilLarge, StencilSmall, wm, input_target_cell, Sorder, Lorder, {0.5,0.5});
 *
 *    printf("The reconstructed value at the given point is %.12f \n", u);
 *
 *    printf("\nThe error measured at this point : %.12f \n\n",abs(u-func(wp->center,{0.0})));
 *
 */
/*
 *    // L2 norm for error analysis
 *    const valarray<double>& gwf = GaussWeightsFace;
 *    const vector< valarray<double> >& gpf = GaussPointsFace;
 *
 *    solution sum;
 *    for (int j=stencilWidth; j<N+stencilWidth; j++){
 *    for (int i=stencilWidth; i<M+stencilWidth; i++){
 *        solution work = 0.0;
 *        vector<valarray<double>> corner;
 *        corner.push_back(mesh[j*(2*stencilWidth+M)+i]);
 *        corner.push_back(mesh[j*(2*stencilWidth+M)+i+1]);
 *        corner.push_back(mesh[(j+1)*(2*stencilWidth+M)+i+1]);
 *        corner.push_back(mesh[(j+1)*(2*stencilWidth+M)+i]);
 *        valarray<double> center;
 *        for (auto & c: corner){
 *            center += c/4.0;
 *        }
 *        point_index target {j,i};
 *        for (int k=0; k<gpf.size(); k++){
 *            valarray<double> mapped = GaussMapPointsFace(gpf[i],corner);
 *            const point tmp = mapped;
 *            double jac = abs(GaussJacobian(gpf[i],corner));
 *            double gw = gwf[i];
 *            work += jac*gw*pow((*func)(mapped,{center[0],center[1]})-WenoPointReconst(StencilLarge, StencilSmall, wm, target, Sorder, Lorder, mapped),2);
 *        }
 *        sum += work;
 *    }}
 *
 *    cout << endl << "L2 norm is calculated as : " << pow(sum,0.5) << endl;
 *
 */

    // Test Numintegral on the edge
/*
 *    vector<point> edge_corner = {{0,0},{1,1}};
 *
 *    double numint_edge = NumIntegralEdge(edge_corner,{0},funcX,funcY);
 *
 *    cout << "Numerical Integral on the given edge " << numint_edge << endl;
 *
 */

    Vec solu;
    solution ** localsolu;
    VecDuplicate(globalu, &solu);
    VecCopy(globalu,solu);

    DMDAVecGetArray(dmu,solu,&localsolu);

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

    // test for 2D Burgers equation
    // Explicit time progression for simplicity
    DrawPressure(dmu, &globalu);   

    double T = 0.5;
    double currentT = 0.0;

    // Spectial case
    double h = 1.0/((double)M*(double)N);

    double dt = 0.1;

    PetscInt       xs,ys,xm,ym;
    ierr = DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

    while(currentT < T){
        currentT += dt;
        const WenoMesh * currentwm = new WenoMesh(M,N,stencilWidth,mesh,localsolu);
        for (int j=ys; j<ys+ym; j++){
        for (int i=xs; i<xs+xm; i++){
            for (int pos = 0; pos<4; pos++){
                point_index target {i-xs+wm->ghost,j-ys+wm->ghost};
                localsolu[j][i] += dt/h * pow(-1,pos+1)*TotalFlux(wm,pos,currentT, StencilLarge,
                                                           StencilSmall, target, Sorder, Lorder,
                                                           funcX, funcY, dfuncX, dfuncY);
            }
        }}
        delete currentwm;
        currentwm = NULL;
    }

    // ============================================================

    DMDAVecRestoreArray(dmu,solu,&localsolu);

    DrawPressure(dmu,&solu);

    DMDAVecRestoreArray(dmu, localu, &lu);

    // Destroy Vectors
    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    VecDestroy(&localu);
    DMDestroy(&dm);

    return 0;
}
