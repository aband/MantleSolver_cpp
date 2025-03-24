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

    return 0.1*point[1];
}

// Location functions
bool interior(const indice& globalCell, 
              const MeshInfo& mi){
   if (globalCell[0] > 0 && globalCell[0] < mi.MPIglobalCellSize[0]-1 &&
       globalCell[1] > 0 && globalCell[1] < mi.MPIglobalCellSize[1]-1){
       return true;
   } else {
       return false;
   }
}

bool edge(const indice& globalCell,
          const MeshInfo& mi){

   // return four edges
   if (// left edge
       (globalCell[0] == 0 && 
        globalCell[1] != 0 && 
        globalCell[1] != mi.MPIglobalCellSize[1]-1) ||
       // right edge
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && 
        globalCell[1] != 0 && 
        globalCell[1] != mi.MPIglobalCellSize[1]-1) ||
       // bottom edge
       (globalCell[1] == 0 && 
        globalCell[0] != 0 && 
        globalCell[0] != mi.MPIglobalCellSize[0]-1) ||
       // top edge
       (globalCell[1] == mi.MPIglobalCellSize[1]-1 && 
        globalCell[0] != 0 && 
        globalCell[0] != mi.MPIglobalCellSize[0]-1) 
      ){
       return true;
   } else {
       return false;
   }

}

bool corner(const indice& globalCell, 
            const MeshInfo& mi){

   // return four corners

   if ((globalCell[0] == 0 && globalCell[1] == 0) ||
       (globalCell[0] == 0 && globalCell[1] == mi.MPIglobalCellSize[1]-1) ||
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && globalCell[1] == 0) ||
       (globalCell[0] == mi.MPIglobalCellSize[0]-1 && globalCell[1] == mi.MPIglobalCellSize[1]-1) 
      ){
       return true;
   } else {
       return false;
   }
}

std::string pfunc(const MeshInfo& mi,
                  const indice& gcell){

    std::string pos;

    if (corner(gcell, mi)){
				pos = "corner";
    } 

    if (edge(gcell, mi)){
pos = "edge";

    }

    if (interior(gcell, mi)){
pos = "interior";

    }
 
    return pos;

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

    int M = 5, N = 10;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);

    double L = 0.05, H = 0.3;
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

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalvec, {(L-h0)/2.0,0.0}, func);
//VecView(globalvec, PETSC_VIEWER_STDOUT_WORLD);
    // Distribute local part to local vectors.
    PetscCall(DMGetLocalVector(dmu, &localvec)); 

    PetscCall(DMGlobalToLocalBegin(dmu, globalvec, INSERT_VALUES, localvec));
    PetscCall(DMGlobalToLocalEnd(dmu, globalvec, INSERT_VALUES, localvec));

    PetscCall(DMDAVecGetArray(dmu, localvec, &locvals));

    // =================================================================

    mluse use = mluse();

    Tensor<double> stencilsol33 = Tensor<double>(2);
    stencilsol33.setSize({3,3});

    Tensor<double> stencilsol22 = Tensor<double>(2);
    stencilsol22.setSize({2,2});

    ml.getsol(stencilsol33, locvals, {0,0}, "(3,3)");

    ml.getsol(stencilsol22, locvals, {0,0}, "(2,2)");

    ml.updatesigma(locvals);

    cout << endl;

    unordered_map<std::string, vector<indice>> interior;
    interior.insert(std::make_pair<std::string, vector<indice>>("(3,3)", { {-1,-1} }));
    interior.insert(std::make_pair<std::string, vector<indice>>("(2,2)", { {-1,-1}, {0,-1}, {0,0}, {-1,0} }));

    use.setmethod("interior", interior);
    use.setbias("interior");

    unordered_map<std::string, vector<indice>> edge;
    edge.insert(std::make_pair<std::string, vector<indice>>
    ("(3,3)", {{0,-1} , {-2,-1}, {-1,0}, {-1,-2}}));
	 edge.insert(std::make_pair<std::string, vector<indice>>
    ("(2,2)", {{-1,-1}, {0,-1} , {0,0} , {-1,0} }));

    use.setmethod("edge", edge);
    use.setbias("edge"); 

    unordered_map<std::string, vector<indice>> corner;
    corner.insert(std::make_pair<std::string, vector<indice>>
    ("(3,3)", {{0,0}, {-2,0}, {0,-2}, {-2,-2}}));
	 corner.insert(std::make_pair<std::string, vector<indice>>
    ("(2,2)", {{-1,-1}, {0,-1} , {0,0} , {-1,0} }));

    use.setmethod("corner", corner);
    use.setbias("corner"); 

    weights wgts_interior;
    weights wgts_edge;
    weights wgts_corner;

    indice target {0,0};
    use.computeWgts("interior", ml, target, h0, wgts_interior);
    use.computeWgts("edge", ml, target, h0, wgts_edge);
    use.computeWgts("corner", ml, target, h0, wgts_corner);

    use.printWgts(wgts_interior);
    use.printWgts(wgts_edge);
    use.printWgts(wgts_corner);

    target = {3,0};
    use.computeWgts("interior", ml, target, h0, wgts_interior);
    use.computeWgts("edge", ml, target, h0, wgts_edge);
    use.computeWgts("corner", ml, target, h0, wgts_corner);

    use.printWgts(wgts_interior);
    use.printWgts(wgts_edge);
    use.printWgts(wgts_corner);

    Tensor<weights> allwgts;

    use.computeWgts(ml, mi, h0, allwgts, pfunc);

    // =================================================================

    DMDAVecRestoreArray(dmu,localvec,&locvals);
    DMRestoreLocalVector(dmu, &localvec); 

    PetscCall(VecDestroy(&globalvec));
 
    PetscCall(VecDestroy(&globalmesh));
    PetscCall(DMDestroy(&dmMesh));
    PetscCall(DMDestroy(&dmu));

    return 1;
}
