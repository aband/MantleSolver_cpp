// Test basis soely
#include <iostream>
#include <petsc.h>
#include "integral.h"
#include "input.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "assemble.h"
#include "util.h"
#include "bndry.h"
#include "solve.h"
#include "error.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

vertex bilinearMap(const vertexSet& corners, const vertex& ref){

    vertex target {0,0};

    target[0] = (1-ref[0])*(1-ref[1])*corners.at(0)[0] + ref[0]*(1-ref[1])*corners.at(1)[0] + ref[0]*ref[1]*corners.at(2)[0] + (1-ref[0])*ref[1]*corners.at(3)[0];
    target[1] = (1-ref[0])*(1-ref[1])*corners.at(0)[1] + ref[0]*(1-ref[1])*corners.at(1)[1] + ref[0]*ref[1]*corners.at(2)[1] + (1-ref[0])*ref[1]*corners.at(3)[1];

    return target;
}

int main(int argc, char ** argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Testing shape functions \n"));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< \n"));

    // =================================================================================

    // Start testing mesh function
    // Initializing problem size with 3X3
    int M = 2, N = 2;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    // Create data management object
    DM    dm;
    Vec   fullmesh;
    const int stencilWidth = 1;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidth, NULL, NULL, &dm);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dm);               CHKERRQ(ierr);
    ierr = DMSetUp(dm);                        CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm, &fullmesh);CHKERRQ(ierr); 

    double L = 1, H = 1;
    double xstart = 0.0, ystart = 0.0;

    ierr = PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL); CHKERRQ(ierr);

    int singleStencilTest = 0;
    double scale = 1;
    ierr = PetscOptionsGetInt(NULL,NULL, "-single", &singleStencilTest, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL, "-scale", &scale, NULL);CHKERRQ(ierr);

    if (singleStencilTest){
        L = L/scale;
        H = H/scale;
        xstart = -L/2.0;
        ystart = -H/2.0;
    }

    MeshParam mp;
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    // Uniform or distorted mesh
    int meshtype=1;
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

    //cout << "Mesh Created. To check full mesh, rerun with -printmesh 1 " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // =================================================================================

    // Contain defined mesh in vector container.
    // and verify it.
    vector< valarray<double> > mesh;
    
    ReadMeshPortion(dm, &fullmesh, mesh);

    //cout << "Converted c array of local mesh into vector container c++ " << endl;
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    // =================================================================================

    DM dmu;

    int cell_ghost = 0;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, M,N, PETSC_DECIDE, PETSC_DECIDE, 1, cell_ghost, NULL, NULL, &dmu);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dmu);               CHKERRQ(ierr);
    ierr = DMSetUp(dmu);                        CHKERRQ(ierr);

    Vec globalu;
    ierr = DMCreateGlobalVector(dmu,&globalu);CHKERRQ(ierr);

    Vec localu; 
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, globalu, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, globalu, INSERT_VALUES, localu);

    // It can be changed later to not be double
    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    // =================================================================================
    // Create MeshInfo object
    MeshInfo mi; 

    // Assign local mesh and local values to mi
    mi.lmesh = mesh;
    mi.localVals = lu;

    AssignValuesMeshInfo(mi,dm,dmu); 

    // =================================================================================
    // Declare classes for shape functions
    basis * basis_   = new basis();
    Hdivmixed * hdiv = new Hdivmixed();
    BRMixed * br     = new BRMixed();

    br->ComputeTotalDOF(mi);
    hdiv->ComputeTotalDOF(mi);

    FILE *fx = fopen("referenceX.dat","w");
    FILE *fy = fopen("referenceY.dat","w");
    FILE *vx = fopen("valx.dat","w");
    FILE *vy = fopen("valy.dat","w");

    int dof=0;
    ierr = PetscOptionsGetInt(NULL,NULL,"-dof",&dof,NULL);CHKERRQ(ierr);

    int seed = 11;
    double h = 1.0/(double) (seed-1) ;

    // Define a mapping from reference [1 0, 0 1] to a quad
    vertex v0 {0.1,-0.2};
    vertex v1 {0.8,0.1};
    vertex v2 {1.2,0.95};
    vertex v3 {-0.05, 1.03};

    vertexSet corners = {v0,v1,v2,v3};

    // reference element
    basis_->GetCorners(mi, {0,0});

    //basis_->GetCorners(corners);
    for (int j=0; j<seed; j++){
    for (int i=0; i<seed; i++){
        vertex sample {i*h, j*h};
//        vertex target ;
        std::array<vertex, 8>  values = hdiv->ComputeHdivmixed(*basis_, sample);
        std::vector<vertex> newvalues = hdiv->EvaluateAll(*basis_, sample);

//        target = bilinearMap(corners, sample);
//        std::array<vertex, 8>  values = hdiv->ComputeHdivmixed(*basis_, target);
//        std::vector<vertex> newvalues = hdiv->EvaluateAll(*basis_, target);
        

        fprintf(fx, "%f ", sample[0]);
        fprintf(fy, "%f ", sample[1]);
        fprintf(vx, "%f ", values[dof][0]);
        fprintf(vy, "%f ", values[dof][1]);
    }}

    fclose(fx);
    fclose(fy);
    fclose(vx);
    fclose(vy);

    // ========================================================================
/*
    std::cout << "boundary dof for BDM element  " << std::endl;
    for (int hdivdof =0; hdivdof<hdiv->getDOF(); hdivdof++){
        std::vector<int> dof = hdiv->GlobalToLocalMapBndry(mi, hdivdof);
        std::cout <<"Current global dof: " << hdivdof << "  ";

        for (const auto& it : dof){
            std::cout << it << "  ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "boundary dof for BR element   " << std::endl;
    for (int brdof =0; brdof<br->getDOF(); brdof++){
        std::vector<int> dof = br->GlobalToLocalMapBndry(mi, brdof);
        std::cout <<"Current global dof: " << brdof << "  ";

        for (const auto& it : dof){
            std::cout << it << "  ";
        }
        std::cout << std::endl;
    }
*/
    // Clear used objects
    DMDAVecRestoreArray(dmu,localu,&lu);
    DMRestoreLocalVector(dmu, &localu); 

    VecDestroy(&fullmesh);
    VecDestroy(&globalu);
    DMDestroy(&dm);
    DMDestroy(&dmu);

    PetscFinalize();

    return 0;
}
