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
/*
 *    if (point[0]<0.50-param[0]){
 *        return sin(point[0])+cos(point[1]);
 *    } else {
 *        return sin(point[0])+cos(point[1])+10;
 *        //return exp(point[0]+point[1]);
 *    }
 *
 */
    return sin(point[0])+cos(point[1]);
    //return point[0]*point[0] + point[1]*point[1];
    //return 0.5;
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

typedef struct{
    DM dm;
    WenoReconst ** wr;
    vector< valarray<double> > mesh;
    int ghost;
    double h;
    vector<int> Lorder;
    vector<int> Sorder;
    index_set StencilLarge;
    vector<index_set> StencilSmall;
    int center_indexl;
    vector<int> center_indexs;
} Ctx;

PetscErrorCode FormFunction(TS ts, PetscReal time, Vec U, Vec F, void * ctx){

    PetscErrorCode ierr;
    Ctx *user = (Ctx*)ctx;
    DM  dm = (DM)user->dm;
    PetscInt M,N,xs,ys,xm,ym,stencilwidth;
    PetscFunctionBeginUser;

    WenoReconst ** wr = user->wr;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    // Get local vector
    Vec localu;
    DMGetLocalVector(dm, &localu);

    DMGlobalToLocalBegin(dm, U, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dm, U, INSERT_VALUES, localu);

    // It can be changed later to not be double
    solution  ** lu;
    DMDAVecGetArray(dm, localu, &lu);

    solution ** f;
    DMDAVecGetArray(dm, F, &f);

    const WenoMesh * currentwm = new WenoMesh(M,N,user->ghost,user->mesh,lu);

    int offset = 4;

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+ym; i++){
        if (j<offset || i<offset || j>N-offset || i>M-offset){
            f[j][i] = lu[j][i];
        } else {
            point_index target {i-xs+user->ghost, j-ys+user->ghost};
            double temp = 0.0;
            for (int pos = 0; pos<4; pos++){
                temp -= 1.0/user->h * TotalFlux(currentwm, pos, time, target, wr, 
                                                funcX, funcY, dfuncX, dfuncY); 
            }
				cout << j << " " << i << " The flux is : " << temp << " h : " << user->h<<  endl;
            f[j][i] = temp;
        }
    }}

    int StencilNum = (M+2)*(N+2);

    for (int s=0; s<StencilNum; s++){
        int j=s/(M+2);
        int i=s%(M+2);
        point_index target {i-1+user->ghost, j-1+user->ghost};
        wr[s] = new WenoReconst(target, currentwm, user->StencilLarge, user->Lorder, user->StencilSmall, user->Sorder, user->center_indexl, user->center_indexs, wr[s]);
        wr[s]->CreateNewWeights2();
    }

    DMDAVecRestoreArray(dm, F, &f);
    DMDAVecRestoreArray(dm, localu, &lu);
    DMRestoreLocalVector(dm, &localu);

    PetscFunctionReturn(0); 
}

/*
 *PetscErrorCode FormJacobian(TS ts, PetscReal time, Vec U, Mat J, Mat Jp, void * ctx){
 *    PetscErrorCode    ierr;
 *    Ctx *user = (Ctx*)ctx;
 *    DM  dm = (DM)user->dm;
 *    PetscInt M,N,xs,ys,xm,ym,stencilwidth;
 *    PetscFunctionBeginUser;
 *
 *    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                                CHKERRQ(ierr);
 *    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);
 *
 *    // Get local vector
 *    Vec localu;
 *    DMGetLocalVector(dm, &localu);
 *
 *    DMGlobalToLocalBegin(dm, U, INSERT_VALUES, localu);
 *    DMGlobalToLocalEnd(dm, U, INSERT_VALUES, localu);
 *
 *    // It can be changed later to not be double
 *    solution  ** lu;
 *    DMDAVecGetArray(dm, localu, &lu);
 *
 *    const WenoMesh * currentwm = new WenoMesh(M,N,user->ghost,user->mesh,lu);
 *
 *    int ns = user->Lorder[0]*user->Lorder[1];
 *
 *    for (int j=ys; j<ys+ym; j++){
 *    for (int i=xs; i<xs+xm; i++){
 *        MatStencil row, col[ns];
 *        row.i=i; row.j = j;
 *        if (j<4 || i<4 || j>N-6 || i>M-6){
 *            col[0].i = i; col[0].j = j; val[0] = 1.0;
 *            nc++;
 *        } else {
 *            point_index target {i-xs+user->ghost, j-ys+user->ghost};
 *
 *            double * df = new double [ns];
 *
 *            for (int pos=0; pos<4; pow++){
 *                double * tempdf = FLuxDerivative(currentwm,pos,time,target,wr,
 *                                                 funcX,funcY,dfuncX,dfuncY);
 *                for (int e=0; e<ns; e++){
 *                    df[e] += -1.0/user->h * tempdf[e];
 *                }
 *            }
 *
 *            // Assign stencil values
 *            for (int e=0; e<ns; e++){
 *                col[e].i = i+user->StencilLarge[e][0];
 *                col[e].j = j+user->StencilLarge[e][1];
 *            }
 *            nc++;
 *
 *        }
 *        MatSetValuesStencil(Jpre, 1, &row, ns, col, df, ADD_VALUES);
 *    }}
 *
 *    MatAssemblyBegin(Jp, MAT_FINAL_ASSEMBLY);
 *    MatAssemblyEnd(Jp, MAT_FINAL_ASSEMBLY);
 *
 *    if (J != Jp){
 *        MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
 *        MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);
 *    }
 *
 *    DMDAVecRestoreArray(dm, localu, &lu);
 *    DMRestoreLocalVector(dm, &localu);
 *
 *    PetscFunctionReturn(0);
 *}
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

    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;


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

    cout << "Mesh Created. To check full mesh, rerun with -printmesh 1 " << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

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

    int center_indexl = 4;

    vector<index_set> StencilSmall;
    StencilSmall.push_back({{-1,-1},{ 0,-1},{ 0, 0},{-1, 0}});
    StencilSmall.push_back({{ 0,-1},{ 1,-1},{ 1, 0},{ 0, 0}});
    StencilSmall.push_back({{ 0, 0},{ 1, 0},{ 1, 1},{ 0, 1}});
    StencilSmall.push_back({{-1, 0},{ 0, 0},{ 0, 1},{-1, 1 }});

    vector<int> center_indexs = {2,3,0,1};

    // Define (5,3) reconstruction
    vector<int> Lorder5 {5,5};
    vector<int> Sorder3 {3,3};

    index_set StencilLarge5 {{-2,-2},{-1,-2},{0 ,-2},{1 ,-2},{2 ,-2},
                             {-2,-1},{-1,-1},{0 ,-1},{1 ,-1},{2 ,-1},
                             {-2, 0},{-1, 0},{0 , 0},{1 , 0},{2 , 0},
                             {-2, 1},{-1, 1},{0 , 1},{1 , 1},{2 , 1},
                             {-2, 2},{-1, 2},{0 , 2},{1 , 2},{2 , 2}};

    int center_indexl5 = 4;

    vector<index_set> StencilSmall3;
    StencilSmall3.push_back({{-2,-2},{-1,-2},{0 ,-2},
                             {-2,-1},{-1,-1},{0 ,-1},
                             {-2, 0},{-1, 0},{0 , 0}});
    StencilSmall3.push_back({{0 ,-2},{1 ,-2},{2 ,-2},
                             {0 ,-1},{1 ,-1},{2 ,-1},
                             {0 , 0},{1 , 0},{2 , 0}});
    StencilSmall3.push_back({{0 , 0},{1 , 0},{2 , 0},
                             {0 , 1},{1 , 1},{2 , 1},
                             {0 , 2},{1 , 2},{2 , 2}});
    StencilSmall3.push_back({{-2, 0},{-1, 0},{0 , 0},
                             {-2, 1},{-1, 1},{0 , 1},
                             {-2, 2},{-1, 2},{0 , 2}});

    vector<int> center_indexs3 = {8,6,0,2};

    // ==========================================================================================================================

    DM dmu;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, 3, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu);               CHKERRQ(ierr);
    ierr = DMSetUp(dmu);                        CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    //SimpleInitialValue(dm,dmu,&fullmesh,&globalu,func);

    // Initialize with oblique data for Burgers equation 
    ObliqueBurgers(dm,dmu,&fullmesh,&globalu,Initial_Condition);
    //SimpleInitialValue(dm,dmu,&fullmesh,&globalu,func);

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

    // Test weno reconstruction
	 point_index target_index_test {N/2+wm->ghost, M/2+wm->ghost};

	 WenoReconst * wr_test;

	 /*
	  *wr_53_test = new WenoReconst(target_index_test,wm,StencilLarge5,Lorder5,StencilSmall3,Sorder3);
	  *wr_53_test->CreateNewWeights();
	  */
	 wr_test = new WenoReconst(target_index_test,wm,StencilLarge,Lorder,StencilSmall,Sorder,center_indexl,center_indexs);
	 wr_test->CreateNewWeights2();

    wr_test->CheckWeights();

	 point target_point_test = {0.5,0.5};

	 printf("\nThe target point is (%.2f, %.2f), the exact value is %.12f \n \n", target_point_test[0], target_point_test[1], func(target_point_test,{0.0}));

	 solution u = wr_test->PointReconstruction(wm, {0.5,0.5});

	 printf("The reconstructed value at the given point is %.12f \n", u);

	 printf("\nThe error measured at this point : %.12f \n\n",abs(u-func(target_point_test,{0.0})));

	 delete wr_test;


    // test for 2D Burgers equation
    // Explicit time progression for simplicity
    DrawPressure(dmu, &globalu);   

    double T = 0.01;
    double currentT = 0.0;

    // Spectial case
    double h = (L*H)/((double)M*(double)N);

    double dt = 0.2*3.0/(double)M;

    PetscInt       xs,ys,xm,ym;
    ierr = DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

    // Create array of weno prepare object
    // First determine how many stencils we need locally
    int StencilNum = (M+2)*(N+2);

    typedef WenoReconst*  wrPtr;

	 wrPtr * wr_23 = new wrPtr[StencilNum];

	 for (int s=0; s<StencilNum; s++){
		  int j=s/(M+2);
		  int i=s%(M+2);
		  point_index target {i-1+wm->ghost, j-1+wm->ghost};
		  wr_23[s] = new WenoReconst(target, wm, StencilLarge, Lorder, StencilSmall, Sorder, center_indexl, center_indexs);
		  wr_23[s]->CreateNewWeights2();
	 }


 //   point_index target {wm->ghost,wm->ghost};

  //  wrPtr wr_23 = new WenoReconst(target,wm,StencilLarge,Lorder,StencilSmall,Sorder,center_indexl,center_indexs);
  //  wr_23->CheckBasisCoeff();

    // Define two types of 4-3 weno reconstruction
/*
 *    wrPtr * wr_53 = new wrPtr[StencilNum];
 *
 *    for (int s=0; s<StencilNum; s++){
 *        int j=s/(M+2);
 *        int i=s%(M+2);
 *        point_index target {i-1+wm->ghost, j-1+wm->ghost};
 *        wr_53[s] = new WenoReconst(target, wm, StencilLarge5, Lorder5, StencilSmall3, Sorder3);
 *        wr_53[s]->CreateNewWeights();
 *        //wr_53[s]->CheckBasisCoeff();
 *        //wr_53[s]->CheckWeights();
 *    }
 *
 */

    cout << "Weno Preparation completed." << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    DMDAVecRestoreArray(dmu,localu,&lu);

/*
 *    solution ** gu;
 *
 *    while(currentT < T){
 *        currentT += dt;
 *        cout << "Current time is : " << currentT << endl;
 *
 *        DMDAVecGetArray(dmu, globalu, &gu);
 *
 *        Vec llu;
 *        DMGetLocalVector(dmu, &llu);
 *        DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, llu);
 *        DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, llu);
 *
 *        DMDAVecGetArray(dmu, llu, &lu);
 *
 *        //Update weno mesh object with new solution propogation through time
 *        const WenoMesh * currentwm = new WenoMesh(M,N,stencilWidth,mesh,lu);
 *        for (int j=ys+3; j<ys+ym-3; j++){
 *        for (int i=xs+3; i<xs+xm-3; i++){
 *            point_index target {i-xs+wm->ghost,j-ys+wm->ghost};
 *            for (int pos = 0; pos<4; pos++){
 *                gu[j][i] -= dt/h * TotalFlux(wm, pos, currentT, target, wr_23,
 *                                             funcX, funcY, dfuncX, dfuncY);
 *            }
 *        }}
 *
 *
 *    for (int s=0; s<StencilNum; s++){
 *        int j=s/(M+2);
 *        int i=s%(M+2);
 *        point_index target {i-1+wm->ghost, j-1+wm->ghost};
 *        //wr_23[s] = new WenoReconst(target, wm, StencilLarge, Lorder, StencilSmall, Sorder);
 *        wr_23[s] = new WenoReconst(target, wm, StencilLarge, Lorder, StencilSmall, Sorder, wr_23[s]);
 *        wr_23[s]->CreateNewWeights();
 *    }
 *
 *        DMDAVecRestoreArray(dmu, llu, &lu);
 *        DMRestoreLocalVector(dmu, &llu);
 *
 *        DMDAVecRestoreArray(dmu, globalu, &gu);
 *
 *    }
 */

	 // Time stepping with TS object
	 TS   ts;
	 SNES snes;
	 Ctx  ctx;

	 // Set up ctx data
	 ctx.wr = wr_23;
	 ctx.Lorder = Lorder;
	 ctx.Sorder = Sorder;
	 ctx.StencilLarge = StencilLarge;
	 ctx.StencilSmall = StencilSmall;
	 ctx.mesh = mesh;
	 ctx.ghost = stencilWidth;
	 ctx.h = h;
	 ctx.dm = dmu;
	 ctx.center_indexl = center_indexl;
	 ctx.center_indexs = center_indexs;

	 TSCreate(PETSC_COMM_WORLD, &ts);
	 TSSetProblemType(ts,TS_NONLINEAR);
	 TSSetType(ts, TSEULER);

	 TSSetMaxTime(ts,T);
	 TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);
	 TSSetDM(ts,dmu);

	 // Customize nonlinear lu[j][i]
	 TSGetSNES(ts,&snes);
	 TSSetTimeStep(ts,dt);
	 TSSetSolution(ts,globalu);

	 TSSetRHSFunction(ts, globalu, FormFunction, &ctx);

	 cout << "Time stepping started." << endl;
	 cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

	 TSSolve(ts,globalu);

	 // ==========================================================================

	 for (int s=0; s<StencilNum; s++){
		  delete wr_23[s];
	 }
	 delete wr_23;

/*
 *    for (int s=0; s<StencilNum; s++){
 *        delete wr_53[s];
 *    }
 *    delete wr_53;
 *
 */
    // ==========================================================================

    DrawPressure(dmu,&globalu);

    hdf5output(dmu,&globalu);

    // Destroy Vectors
    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    DMDestroy(&dm);
//    TSDestroy(&ts);


    return 0;
}
