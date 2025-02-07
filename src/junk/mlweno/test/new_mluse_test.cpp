// A new test implementing new MLWENO weighting and field 
// And testing reconstruction near boundary 
// Using default stencil management in driver file

#include "driver.h"

int main(int argc, char ** argv){

    PetscMPIInt size, rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Input mesh parameter =========================================================
    int M=3, N=3;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL));

    double L = 0.2, H = 1;
    double xstart = -0.1, ystart = -1.0001;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    int stencilWidthMesh = 5;
    int stencilWidthU = 3;

    int physicsScale = 1;
    PetscCall(PetscOptionsGetInt(NULL,NULL, "-scale", &physicsScale, NULL));

    int meshType = 0; 
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshType,NULL));

    int maxIter = 1; 
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-maxIter", &maxIter, NULL));        
    double tolUzawa = 10e-7; 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tol", &tolUzawa, NULL)); 

    double Tmax = 0.009; // Stop at the first step 
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tmax", &Tmax, NULL)); 

    int targeti=1, targetj=1;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-i",&targeti,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-j",&targetj,NULL));

    // ==============================================================================

    Driver * driver = new Driver();

    driver->CreatePhase();

    driver->CreateMesh(M, N, L, H, xstart, ystart, 
                       stencilWidthMesh, stencilWidthU,
                       physicsScale, meshType);

    driver->InitTransport(InitHD, InitCD);

//    PetscCall(VecView(driver->globalHD, PETSC_VIEWER_STDOUT_WORLD));
//    PetscCall(VecView(driver->globalCD, PETSC_VIEWER_STDOUT_WORLD));

    driver->PrepareDefaultTransport();

    Vec localCD, localHD;

    PetscCall(DMGetLocalVector(driver->dmu, &localCD));
    PetscCall(DMGetLocalVector(driver->dmu, &localHD));

    PetscCall(DMGlobalToLocalBegin(driver->dmu, driver->globalCD, INSERT_VALUES, localCD));
    PetscCall(DMGlobalToLocalEnd(driver->dmu, driver->globalCD, INSERT_VALUES, localCD));
    PetscCall(DMGlobalToLocalBegin(driver->dmu, driver->globalHD, INSERT_VALUES, localHD));
    PetscCall(DMGlobalToLocalEnd(driver->dmu, driver->globalHD, INSERT_VALUES, localHD));

    double ** lc;
    double ** lh;

    DMDAVecGetArray(driver->dmu, localCD, &lc);
    DMDAVecGetArray(driver->dmu, localHD, &lh);

    //! Hook double pointer to meshInfo object for the calculation of smoothness indicator 
	 //! and non linear weights later

    driver->mi.localCD = lc;
    driver->mi.localHD = lh;

    driver->UpdateSmoothnessIndicator();
    driver->UpdateNonlinearWgts();

    const valarray<double>& gpe = GaussPointsEdge;

    // Compute flux on one given edge
    indice target {targeti,targetj};

    // Perform reconstruction in a controlled manner 
    cout << "Let's see what actually happens at four edges of cell ("<<targeti << 
				", " << targetj << ") : " << endl;

    vertexSet corners = extractCorners(driver->mi, target);

    std::vector<vertex> gauss_p;
    gauss_p.resize(gpe.size());

    // Loop through four edges
    for (int e=0; e<4; e++){
        std::vector<vertex> edge {corners.at((e+3)%4), corners.at(e)};
        for (int g=0; g<gpe.size(); g++){
            vertex rp = GaussMapPointsEdge({gpe[g]}, edge);
            double val = driver->mluseAdv_->Evaluate(rp, target, driver->mi, location(driver->mi, target), "HD");
            cout <<val<< " ";
        }
        cout << endl;
    }

    // Restore HD and CD
    // (Update ghost region)
    // global vector not touched 
    PetscCall(DMDAVecRestoreArray(driver->dmu, localCD, &lc));
    PetscCall(DMDAVecRestoreArray(driver->dmu, localHD, &lh));
    PetscCall(DMRestoreLocalVector(driver->dmu, &localCD));
    PetscCall(DMRestoreLocalVector(driver->dmu, &localHD));

    return 1;
}
