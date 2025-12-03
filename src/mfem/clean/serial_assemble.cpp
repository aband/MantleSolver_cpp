#include "serial_solver.h"

static int ExtractCellPorosity(const vector<double>& edgeporo, 
                               const vector<double>& cellporo,
                               const vector<double>& averporo,
                               int i, int j, int M, int N, 
                               poroSet& poro){

    poro.aveporo = averporo.at(j*M+i);

    // Extract porosity on cell gauss points
    poro.cellporo.clear();
    poro.cellporo.resize(9);

    for (int g=0; g<9; g++){
        poro.cellporo.at(g) = cellporo.at((j*M+i)*9 + g);
    }

    // Extract porosity on edge gauss points
    poro.edgeporo.clear();
    poro.edgeporo.resize(12);

    // Create index set for four edges
    vector<int> indexSet {(j*(M_+1) + i)*3  , (j*M_     + i)*3, 
                          (j*(M_+1) + i+1)*3, ((j+1)*M_ + i)*3}; 

    for (int e=0; e<4; e++){
    for (int g=0; g<3; g++){
        poro.edgeporo.at(e*3+g) = edgeporo.at(indexSet.at(e)+g);
    }}

    return 1;
}

inline bool elemOnBndry(const MeshInfo& mi,
                        const indice& global){
     
    if (global[0] == 0 ||
        global[1] == 0 ||
        global[0] == mi.MPIglobalCellSize[0] - 1 ||
        global[1] == mi.MPIglobalCellSize[1] - 1){

        return true;

    } else {

        return false;

    }
}

int DarcyStokes::PrepareReducedSys(int Adnz, int Aonz, int Bdnz, int Bonz){

    // Create space holding reduced linear matrices
    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, reducedDOF, 
                                             Adnz, NULL, Aonz, NULL, &redsys->M));
    PetscCall(MatSetUp(redsys->M));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, bndrySize,
                                             Adnz, NULL, Aonz, NULL, &redsys->Kg));  
    PetscCall(MatSetUp(redsys->Kg));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, totalElem,
                                             Bdnz, NULL, Bonz, NULL, &redsys->B));   
    PetscCall(MatSetUp(redsys->B));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             bndrySize, totalElem,
                                             Bdnz, NULL, Bonz, NULL, &redsys->Bg));  
    PetscCall(MatSetUp(redsys->Bg));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             totalElem, totalElem,
                                             1, NULL, 0, NULL, &redsys->C));   
    PetscCall(MatSetUp(redsys->C));

    // Create Corresponding vectors
    // Get ownership first
    // g vector has the size of bndrySize
    // source vector has the size of reducedDOF 
    // neum vector has the sie of reducedDOF the same as source vector
    int m, n; 
    PetscCall(MatGetOwnershipRange(redsys->Bg, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys->g));
    PetscCall(VecSetUp(redsys->g));

    PetscCall(MatGetOwnershipRange(redsys->B, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys->source));
    PetscCall(VecSetUp(redsys->source));

    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys->neum));
    PetscCall(VecSetUp(redsys->neum));

    return 1;
}

int DarcyStokes::Assemble(const MeshInfo& mi,
                          const vector<double>& edgeporo,
                          const vector<double>& cellporo,
                          const vector<double>& averporo,
                          double theta){

    PetscMPIInt size, rank;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscFunctionBeginUser;

    // Get gauss points first
//    const valarray<double>& gwe = GaussWeightsEdge;
//    const valarray<double>& gpe = GaussPointsEdge;
//    const valarray<double>& gwf = GaussWeightsFace;
//    const vector<vertex>&   gpf = GaussPointsFace;

    // Calculate dofs 
    totalElem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    reducedDOFStokes = br_->getDOF() - bndryDOFStokes_;

    reducedDOFDarcy = hdiv_->getDOF() - bndryDOFDarcy_;

    // Initialize reduced linear system
    PrepareReducedSys(reducedStokes_, reducedDOFStokes, bndryDOFStokes_, 
                      totalElem, 30, 30, 4, 4);
    PrepareReducedSys(reducedDarcy_, reducedDOFDarcy, bndryDOFDarcy_, 
                      totalElem, 14, 14, 2, 2);

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           totalElem, totalElem, 
                           1, NULL, 0, NULL, &K));
    PetscCall(MatSetUp(K));

    // =====================================================================
    LocMat * locmatS = new LocMat;
    LocMat * locmatD = new LocMat;

    double k = 0.0;

    poroSet locporo; 

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};
        basis_.GetCorners(mi, global);

        ExtractCellPorosity(edgeporo, cellporo, averporo, 
                            i, j, mi.MPIglobalCellSize[0], mi.MPIglobalCellSize[1], 
                            locporo);        

        AssignLocMatStokes(mi, locmatS, theta, locporo);
        AssignLocMatDarcy(mi, locmatD, theta, locporo);
        AssignLocMatCouple(mi, k, theta, locporo);

        if (elemOnBndry(mi, global)){

            AssignLocRedSysBndry(redsysStokes, locmatS, refArrayStokes, mi, 
                                 bndryEssenStokesAll, );
            AssignLocRedSysBndry();

        } else {

            AssignLocRedSys();
            AssignLocRedSys();
        }

        // Assign coupling K matrix and two C matrices
        // const pressure space not affected by boundary dofs
        PetscCall(MatSetValue(K,nElem,nElem,k,ADD_VALUES));
        PetscCall(MatSetValue(reducedStokes_->C, nElem, nElem, locmatS->C,ADD_VALUES));
        PetscCall(MatSetValue(reducedDarcy_->C, nElem, nElem, locmatD->C, ADD_VALUES));

    }}

    return 1;
}
