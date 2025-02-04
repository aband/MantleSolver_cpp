#include "stencilpolynomial.h"
#include "reconstruction.h"
#include "mluse.h"
#include "petsc.h"
#include "input.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

double func(const vertex& point,
            const vector<double>& param){

    if (point[0] < param[0]) {

    return sin(point[0])*cos(point[1]);

    } else {

    return sin(point[0])*cos(point[1]) + 0.0;

    }
}

int main(int argc, char ** argv){

    // Integrated test with petsc and mesh functions

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    int M = 5, N = 5;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

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

    double h0 = sqrt((L*H)/(double)(M*N));

    vertex test {0.5,0.5};

    test = (test-h0 /3)/dscale;

    multilevel ml = multilevel();

    ml.addLevel("(3,3)", {3,3}, mi);
    ml.addLevel("(2,2)", {2,2}, mi);

    // =================================================================
    Vec globalvec, localvec;
    double ** locvals;

    PetscCall(DMCreateGlobalVector(dmu, &globalvec));

//    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalvec, {(L-h0)/2.0,0.0}, func);

    const PetscScalar y[25] = { 0.271849 , 0.994933, 0.408574, 0.00365381, -3.21759e-20, 0.272165, 0.993494, 0.410133, 0.00322378,-3.34346e-18, 0.272162, 0.993547, 0.410113, 0.00320021, -3.66537e-18, 0.271801, 0.996462 ,0.407558 ,0.00317677 ,-3.1735e-18 ,0.271842 ,0.994967 ,0.408554 ,0.00365379 ,-2.84709e-20 };
    const PetscInt id[25] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    VecSetValues(globalvec, 25, id, y, INSERT_VALUES);

   // VecView(globalvec, PETSC_VIEWER_STDOUT_WORLD);

    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvec)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalvec, INSERT_VALUES, localvec));
    PetscCall(DMGlobalToLocalEnd(dmu, globalvec, INSERT_VALUES, localvec));

    PetscCall(DMDAVecGetArray(dmu, localvec, &locvals));

    // =================================================================

    //Updating sigma , dsigma , scaled sigma and d scaled sigma all at once
    ml.updateall(locvals, h0,1,1e-4, mi);

/*
    ml.printsigma("(3,3)");
    ml.printsigma("(2,2)");

    ml.printscaledsigma("(3,3)");
    ml.printscaledsigma("(2,2)");

    ml.printdsigma("(3,3)");
    ml.printdsigma("(2,2)");

    ml.printdscaledsigma("(3,3)");
    ml.printdscaledsigma("(2,2)");
*/

    mluse use = mluse();

    unordered_map<std::string, vector<indice>> method;
    method.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-1,-1} }));
    method.insert(std::make_pair<std::string, vector<indice>>("(2,2)", { {-1,-1}, {0,-1}, {0,0}, {-1,0} }));

    use.setmethod("all",method);
    use.setbias("all");

    Tensor<weights> allwgts;

    use.computeWgts(ml, mi, allwgts);

    indice target {M/2,N/2};

    cout << "Reconstructed value : " << use.eval(test, ml, "all", 
         allwgts({target[0], target[1]}), target, locvals) << endl 
         << "Function value : " << func(test, {(L-h0)/2,0.0})<< endl;

    derivative testder;
//    use.der(test, ml, "all", allwgts({target[0], target[1]}), {target[0],target[1]}, locvals,
//             mi, testder); 

    derivative sum;

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
    use.der(test, ml, "all", allwgts({i, j}), {i,j}, locvals,
             mi, testder); 

    unordered_map_arithmetic(sum, testder, std::plus<double>()); 
    }}

//    use.der(test, ml, "all", allwgts({1, 2}), {1,2}, locvals,
//             mi, testder); 

    unordered_map_print(sum);
    // =================================================================

    DMDAVecRestoreArray(dmu,localvec,&locvals);
    DMRestoreLocalVector(dmu, &localvec); 

    PetscCall(VecDestroy(&globalvec));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    return 1;
}

