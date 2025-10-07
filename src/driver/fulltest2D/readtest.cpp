#include "driver.h"
#include "print.h"
#include "read.h"

int main(int argc, char **argv){

    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Input mesh parameter =========================================================
    int M=5, N=20;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL));

    double L = 0.1, H = 0.4;
    double xstart = -0.5*L, ystart = -1.0001*H;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    int stencilWidthMesh = 5;
    int stencilWidthU = 3;

    int physicsScale = 0;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-scale", &physicsScale, NULL));

    int meshType = 0; 
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshType,NULL));

    int maxIter = 20; 
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL));        
    double tolUzawa = 10e-17; 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL)); 

    double Tmax = 20; // Stop at the first step 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tmax", &Tmax, NULL)); 

    double dt = 1;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));

    int showPhase = 0;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-showphase", &showPhase, NULL)); 

    int withUnit = 0;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-unit", &withUnit, NULL));

    int interval = 1;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-interval", &interval, NULL));

    // ==============================================================================

    Driver * driver = new Driver();

    driver->withUnit = withUnit;
    driver->CreatePhase();
    driver->ShowPhase();
    driver->CreateMesh(M, N, L, H, xstart, ystart, 
                       stencilWidthMesh, stencilWidthU,
                       physicsScale, meshType);

    std::vector<double> restartHD; restartHD.resize(M*N);
    std::vector<double> restartCD; restartCD.resize(M*N);

    ReadValues("restartHD.dat", restartHD);
    ReadValues("restartCD.dat", restartCD);

    PetscCall(DMCreateGlobalVector(driver->dmu, &driver->globalCD));
    PetscCall(DMCreateGlobalVector(driver->dmu, &driver->globalHD));

    SimpleInitialValue(driver->dmu, &driver->globalCD, restartCD);
    SimpleInitialValue(driver->dmu, &driver->globalHD, restartHD);

    multilevel ml = multilevel();
    // 1D test ================================================================================================
    ml.addLevel("(1,3)", {1,3}, driver->mi);
    ml.addLevel("(1,2)", {1,2}, driver->mi);

    //ml.addLevel("(3,3)", {3,3}, driver->mi);
    //ml.addLevel("(2,2)", {2,2}, driver->mi);

    // Area scale
    double h0 = sqrt((L*H)/
         (double)(driver->mi.MPIglobalCellSize[0]*driver->mi.MPIglobalCellSize[1]));

    mluse testuse = mluse();

    unordered_map<std::string, vector<indice>> method;
    method.insert(std::make_pair<std::string, vector<indice>>("(1,3)", { {0,-1} }));
    method.insert(std::make_pair<std::string, vector<indice>>("(1,2)", { {0,-1} , {0,0} }));

    testuse.setmethod("test", method);
    testuse.setbias("test");

    // 2D test ================================================================================================
/*
    unordered_map<std::string, vector<indice>> interior;
    interior.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-1,-1} }));
    interior.insert(std::make_pair<std::string, vector<indice>>("(2,2)", { {-1,-1}, {0,-1}, {0,0}, {-1,0} }));

    testuse.setmethod("interior", interior);
    testuse.setbias("interior");

    unordered_map<std::string, vector<indice>> edge;
    edge.insert(std::make_pair<std::string, vector<indice>>
    ("(3,3)", {{0,-1} , {-2,-1}, {-1,0}, {-1,-2}}));
	 edge.insert(std::make_pair<std::string, vector<indice>>
    ("(2,2)", {{-1,-1}, {0,-1} , {0,0} , {-1,0} }));

    testuse.setmethod("edge", edge);
    testuse.setbias("edge"); 

    unordered_map<std::string, vector<indice>> corner;
    corner.insert(std::make_pair<std::string, vector<indice>>
    ("(3,3)", {{0,0}, {-2,0}, {0,-2}, {-2,-2}}));
	 corner.insert(std::make_pair<std::string, vector<indice>>
    ("(2,2)", {{-1,-1}, {0,-1} , {0,0} , {-1,0} }));

    testuse.setmethod("corner", corner);
    testuse.setbias("corner"); 
*/


    // Get local data arrays and vectors
    Vec localHD, localCD;
    double ** lHD;
    double ** lCD;

    PetscCall(DMGetLocalVector(driver->dmu, &localHD)); 

    PetscCall(DMGlobalToLocalBegin(driver->dmu, driver->globalHD, INSERT_VALUES, localHD));
    PetscCall(DMGlobalToLocalEnd(driver->dmu, driver->globalHD, INSERT_VALUES, localHD));

    PetscCall(DMDAVecGetArray(driver->dmu, localHD, &lHD));

    PetscCall(DMGetLocalVector(driver->dmu, &localCD)); 

    PetscCall(DMGlobalToLocalBegin(driver->dmu, driver->globalCD, INSERT_VALUES, localCD));
    PetscCall(DMGlobalToLocalEnd(driver->dmu, driver->globalCD, INSERT_VALUES, localCD));

    PetscCall(DMDAVecGetArray(driver->dmu, localCD, &lCD));

    // ==============================================================

    ml.updatesigma(lHD);
    Tensor<weights> allwgtsHD;
    testuse.computeWgts(ml, driver->mi, h0, allwgtsHD, location);

    ml.updatesigma(lCD);
    Tensor<weights> allwgtsCD;
    testuse.computeWgts(ml, driver->mi, h0, allwgtsCD, location);

    for (int j=0; j<N; j++){
        for (int i=0; i<M; i++){
            //testuse.printWgts(allwgtsHD({i,j}));
            testuse.printWgts(allwgtsCD({i,j}));
        }cout << endl;
    }

    // ==============================================================

    DMDAVecRestoreArray(driver->dmu,localHD,&lHD);
    DMRestoreLocalVector(driver->dmu, &localHD); 
    DMDAVecRestoreArray(driver->dmu,localCD,&lCD);
    DMRestoreLocalVector(driver->dmu, &localCD); 

    return 1;
}
