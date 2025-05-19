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

    FILE *fx1 = fopen("gridX1.dat","w");
    FILE *fy1 = fopen("gridY1.dat","w");
    FILE *vx1 = fopen("valx1.dat","w");
    FILE *vy1 = fopen("valy1.dat","w");
    FILE *ux1 = fopen("ualx1.dat","w");
    FILE *uy1 = fopen("ualy1.dat","w");
    FILE *sc1 = fopen("sc1.dat", "w");
    FILE *scc1 = fopen("scc1.dat", "w");

    FILE *fx2 = fopen("gridX2.dat","w");
    FILE *fy2 = fopen("gridY2.dat","w");
    FILE *vx2 = fopen("valx2.dat","w");
    FILE *vy2 = fopen("valy2.dat","w");
    FILE *ux2 = fopen("ualx2.dat","w");
    FILE *uy2 = fopen("ualy2.dat","w");
    FILE *sc2 = fopen("sc2.dat", "w");
    FILE *scc2 = fopen("scc2.dat", "w");


    int dof=0;
    ierr = PetscOptionsGetInt(NULL,NULL,"-dof",&dof,NULL);CHKERRQ(ierr);

    int seed = 101;
    double h = 1.0/(double) (seed-1) ;

    // Define a mapping from reference [1 0, 0 1] to a quad
    vertex v0 {0.1,-0.2};
    vertex v1 {0.8,0.1};
    vertex v2 {1.2,0.95};
    vertex v3 {-0.05, 1.03};

    vertexSet corners1 = {v0,v1,v2,v3};

    vertex vv0 {0.8,0.1};
    vertex vv1 {1.65,0.02};
    vertex vv2 {1.7,1.0};
    vertex vv3 {1.2,0.95};

    vertexSet corners2 = {vv0,vv1,vv2,vv3};

    // reference element
    //basis_->GetCorners(mi, {0,0});
        vertexSet edge {v1,v2};
        vertex unitNormal = UnitNormal(edge, length(edge));

    //basis_->GetCorners(corners);
    for (int j=0; j<seed; j++){
    for (int i=0; i<seed; i++){
        vertex sample {i*h, j*h};
//        vertex target ;
//        std::array<vertex, 8>  values = hdiv->ComputeHdivmixed(*basis_, sample);
//        std::vector<vertex> newvalues = hdiv->EvaluateAll(*basis_, sample);

        basis_->GetCorners(corners1);
        vertex target1 = bilinearMap(corners1, sample);
        std::array<vertex, 8>  values1 = hdiv->ComputeHdivmixed(*basis_, target1);
        //std::vector<vertex> newvalues1 = hdiv->EvaluateAll(*basis_, target);
	     std::array<vertex, 12> brvalues1 = br->ComputeBRmixed(*basis_, target1);
     
        basis_->GetCorners(corners2);
        vertex target2 = bilinearMap(corners2, sample);
        std::array<vertex, 8>  values2 = hdiv->ComputeHdivmixed(*basis_, target2);
 	     std::array<vertex, 12> brvalues2 = br->ComputeBRmixed(*basis_, target2);
 
        fprintf(fx1, "%f ", target1[0]);
        fprintf(fy1, "%f ", target1[1]);
        fprintf(fx2, "%f ", target2[0]);
        fprintf(fy2, "%f ", target2[1]);
/*
        fprintf(vx1, "%f ", values1[2][0]);
        fprintf(vy1, "%f ", values1[2][1]);
        fprintf(vx2, "%f ", values2[0][0]);
        fprintf(vy2, "%f ", values2[0][1]);

        fprintf(ux1, "%f ", values1[2+4][0]);
        fprintf(uy1, "%f ", values1[2+4][1]);
        fprintf(ux2, "%f ", values2[0+4][0]);
        fprintf(uy2, "%f ", values2[0+4][1]);
*/
        fprintf(vx1, "%f ", brvalues1[8+2][0]);
        fprintf(vy1, "%f ", brvalues1[8+2][1]);
        fprintf(vx2, "%f ", brvalues2[8+0][0]);
        fprintf(vy2, "%f ", brvalues2[8+0][1]);

        fprintf(ux1, "%f ", brvalues1[2+4][0]);
        fprintf(uy1, "%f ", brvalues1[2+4][1]);
        fprintf(ux2, "%f ", brvalues2[0+4][0]);
        fprintf(uy2, "%f ", brvalues2[0+4][1]);
/*
        fprintf(sc1, "%f ", values1[2][0]*unitNormal[0] + values1[2][1]*unitNormal[1]);
        fprintf(sc2, "%f ", values2[0][0]*unitNormal[0] + values2[0][1]*unitNormal[1]);
		  */
        fprintf(sc1, "%f ", brvalues1[8+2][0]*unitNormal[0] + brvalues1[8+2][1]*unitNormal[1]-1);
        fprintf(sc2, "%f ", brvalues2[8+0][0]*unitNormal[0] + brvalues2[8+0][1]*unitNormal[1]-1);

        fprintf(scc1, "%f ", values1[2+4][0]*unitNormal[0] + values1[2+4][1]*unitNormal[1]);
        fprintf(scc2, "%f ", values2[0+4][0]*unitNormal[0] + values2[0+4][1]*unitNormal[1]);

    }}

    fclose(fx1);
    fclose(fy1);
    fclose(vx1);
    fclose(vy1);
    fclose(ux1);
    fclose(uy1);
    fclose(fx2);
    fclose(fy2);
    fclose(vx2);
    fclose(vy2);
    fclose(ux2);
    fclose(uy2);
    fclose(sc1);
	 fclose(sc2);
    fclose(scc1);
	 fclose(scc2);

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
