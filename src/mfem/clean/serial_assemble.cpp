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
    vector<int> indexSet {(j*(M+1) + i)*3  , (j*M     + i)*3, 
                          (j*(M+1) + i+1)*3, ((j+1)*M + i)*3}; 

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

int DarcyStokes::PrepareReducedSys(ReducedSys& redsys,
                                   int reducedDOF, int bndrySize, 
                                   int Adnz, int Aonz, int Bdnz, int Bonz){

    // Create space holding reduced linear matrices
    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, reducedDOF, 
                                             Adnz, NULL, Aonz, NULL, &redsys.M));
    PetscCall(MatSetUp(redsys.M));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, bndrySize,
                                             Adnz, NULL, Aonz, NULL, &redsys.Kg));  
    PetscCall(MatSetUp(redsys.Kg));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             reducedDOF, totalElem,
                                             Bdnz, NULL, Bonz, NULL, &redsys.B));   
    PetscCall(MatSetUp(redsys.B));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             bndrySize, totalElem,
                                             Bdnz, NULL, Bonz, NULL, &redsys.Bg));  
    PetscCall(MatSetUp(redsys.Bg));

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                                             totalElem, totalElem,
                                             1, NULL, 0, NULL, &redsys.C));   
    PetscCall(MatSetUp(redsys.C));

    // Create Corresponding vectors
    // Get ownership first
    // g vector has the size of bndrySize
    // source vector has the size of reducedDOF 
    // neum vector has the sie of reducedDOF the same as source vector
    int m, n; 
    PetscCall(MatGetOwnershipRange(redsys.Bg, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys.g));
    PetscCall(VecSetUp(redsys.g));

    PetscCall(MatGetOwnershipRange(redsys.B, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys.source));
    PetscCall(VecSetUp(redsys.source));

    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys.neum));
    PetscCall(VecSetUp(redsys.neum));

    return 1;
}

// Used for interior elements (no need to identify boundary dofs)
// Set multiple values at the same time
// need const int id array
// Stokes and Darcy part are separated

int DarcyStokes::AssignLocRedSysDarcy(LocMat& loc,
                                      int * ref,
                                      const MeshInfo& mi,
                                      const indice& global){

    const int idxm = FlatIndic(mi, global);

    std::array<int, 8> locDof = hdiv_.LocalToGlobal(mi, global);

    const int IS[8] = {ref[locDof[0]], ref[locDof[1]], 
                       ref[locDof[2]], ref[locDof[3]],
                       ref[locDof[4]], ref[locDof[5]], 
                       ref[locDof[6]], ref[locDof[7]]};

    const double locB[8] = {loc.B.at(0), loc.B.at(1),
                            loc.B.at(2), loc.B.at(3),
                            loc.B.at(4), loc.B.at(5),
                            loc.B.at(6), loc.B.at(7)};

    // Only fill B and M matrix
    PetscCall(MatSetValues(reducedDarcy_.B, 8, IS, 1, &idxm, locB, ADD_VALUES));

    for (unsigned int l=0; l<8; l++){
        const double locA[8] = {loc.A.at(0+l*8), loc.A.at(1+l*8),
                                loc.A.at(2+l*8), loc.A.at(3+l*8),
                                loc.A.at(4+l*8), loc.A.at(5+l*8),
                                loc.A.at(6+l*8), loc.A.at(7+l*8)};
        const int Aidxm = IS[l];
        PetscCall(MatSetValues(reducedDarcy_.M,1,&Aidxm,8,IS,locA,ADD_VALUES));
    }

    // Fill Source vector
    const double locf[8] = {loc.f.at(0), loc.f.at(1),
                            loc.f.at(2), loc.f.at(3),
                            loc.f.at(4), loc.f.at(5),
                            loc.f.at(6), loc.f.at(7)};

    PetscCall(VecSetValues(reducedDarcy_.source, 8, IS, locf, ADD_VALUES));

    return 0;
}

int DarcyStokes::AssignLocRedSysStokes(LocMat& loc,
                                       int * ref,
                                       const MeshInfo& mi,
                                       const indice& global){

    const int idxm = FlatIndic(mi, global);

    std::array<int, 12> locDof = br_.LocalToGlobal(mi, global);

    const int IS[12] = {ref[locDof[0]], ref[locDof[1]], 
                        ref[locDof[2]], ref[locDof[3]],
                        ref[locDof[4]], ref[locDof[5]], 
                        ref[locDof[6]], ref[locDof[7]],
                        ref[locDof[8]], ref[locDof[9]],
                        ref[locDof[10]], ref[locDof[11]]};

    const double locB[12] = {loc.B.at(0), loc.B.at(1),
                             loc.B.at(2), loc.B.at(3),
                             loc.B.at(4), loc.B.at(5),
                             loc.B.at(6), loc.B.at(7),
                             loc.B.at(8), loc.B.at(9),
                             loc.B.at(10), loc.B.at(11)};

    // Only fill B and M matrix
    PetscCall(MatSetValues(reducedStokes_.B, 12, IS, 1, &idxm, locB, ADD_VALUES));

    for (unsigned int l=0; l<12; l++){
        const double locA[12] = {loc.A.at(0+l*12), loc.A.at(1+l*12),
                                 loc.A.at(2+l*12), loc.A.at(3+l*12),
                                 loc.A.at(4+l*12), loc.A.at(5+l*12),
                                 loc.A.at(6+l*12), loc.A.at(7+l*12),
                                 loc.A.at(8+l*12), loc.A.at(9+l*12),
                                 loc.A.at(10+l*12), loc.A.at(11+l*12)};
        const int Aidxm = IS[l];
        PetscCall(MatSetValues(reducedStokes_.M,1,&Aidxm,12,IS,locA,ADD_VALUES));
    }

    // Only contribute to source vector not g vector

    const double locf[12] = {loc.f.at(0), loc.f.at(1),
                             loc.f.at(2), loc.f.at(3),
                             loc.f.at(4), loc.f.at(5),
                             loc.f.at(6), loc.f.at(7),
                             loc.f.at(8), loc.f.at(9),
                             loc.f.at(10), loc.f.at(11)};

    PetscCall(VecSetValues(reducedStokes_.source, 12, IS, locf, ADD_VALUES));

    return 0;
}

inline int AssembleReducedSys(ReducedSys& redsys){

    PetscCall(MatAssemblyBegin(redsys.M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys.M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(redsys.Kg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys.Kg, MAT_FINAL_ASSEMBLY));

    PetscCall(MatAssemblyBegin(redsys.B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys.B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(redsys.Bg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys.Bg, MAT_FINAL_ASSEMBLY));

    PetscCall(MatAssemblyBegin(redsys.C, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys.C, MAT_FINAL_ASSEMBLY));

    PetscCall(VecAssemblyBegin(redsys.g));
    PetscCall(VecAssemblyEnd(redsys.g));

    PetscCall(VecAssemblyBegin(redsys.source));
    PetscCall(VecAssemblyEnd(redsys.source));

    PetscCall(VecAssemblyBegin(redsys.neum));
    PetscCall(VecAssemblyEnd(redsys.neum));

    return 1;
}

int DarcyStokes::Assemble(const MeshInfo& mi,
                          const vector<double>& edgeporo,
                          const vector<double>& cellporo,
                          const vector<double>& averporo,
                          double theta,
                          PhysProperty * pp){

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

    reducedDOFStokes = br_.getDOF() - bndryDOFStokesEssen_;

    reducedDOFDarcy = hdiv_.getDOF() - bndryDOFDarcyEssen_;

    // Initialize reduced linear system
    PrepareReducedSys(reducedStokes_, reducedDOFStokes, bndryDOFStokesEssen_, 
                      30, 30, 4, 4);
    PrepareReducedSys(reducedDarcy_, reducedDOFDarcy, bndryDOFDarcyEssen_, 
                      14, 14, 2, 2);

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           totalElem, totalElem, 
                           1, NULL, 0, NULL, &K));
    PetscCall(MatSetUp(K));

    // =====================================================================
    LocMat locmatS;
    LocMat locmatD;

    double k = 0.0;

    poroSet locporo; 

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice global {i,j};
        int nElem = FlatIndic(mi, global); 
        basis_.GetCorners(mi, global);

        ExtractCellPorosity(edgeporo, cellporo, averporo, 
                            i, j, mi.MPIglobalCellSize[0], mi.MPIglobalCellSize[1], 
                            locporo);        

        AssignLocMatStokes(mi, locmatS, theta, locporo, pp);
        AssignLocMatDarcy(mi, locmatD, theta, locporo, pp);
        AssignLocMatCouple(mi, k, theta, locporo, pp);

        if (elemOnBndry(mi, global)){

            AssignLocRedSys(reducedStokes_, locmatS, refArrayStokesEssen_, refArrayStokesNatur_,  
                            mi, bndryStokesAll, global, br_, {0.0});
            AssignLocRedSys(reducedDarcy_, locmatD, refArrayDarcyEssen_, refArrayDarcyNatur_,
                            mi, bndryDarcyAll, global, hdiv_, {0.0});

        } else {

            AssignLocRedSysStokes(locmatS, refArrayStokesEssen_, mi, global);
            AssignLocRedSysDarcy( locmatD, refArrayDarcyEssen_ , mi, global);
        }

        // Assign coupling K matrix and two C matrices
        // const pressure space not affected by boundary dofs
        PetscCall(MatSetValue(K,nElem,nElem,k,ADD_VALUES));
        PetscCall(MatSetValue(reducedStokes_.C, nElem, nElem, locmatS.C,ADD_VALUES));
        PetscCall(MatSetValue(reducedDarcy_.C, nElem, nElem, locmatD.C, ADD_VALUES));

    }}

    AssembleReducedSys(reducedStokes_);
    AssembleReducedSys(reducedDarcy_);

    PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

    return 1;
}

static int RearrangeLinearSys(ReducedSys& redsys, const int& nelem){

    // Right hand side F
    PetscCall(VecDuplicate(redsys.source, &redsys.F));
    PetscCall(MatMult(redsys.Kg, redsys.g, redsys.F));

    PetscCall(VecAYPX(redsys.F, -1, redsys.source));

    // Right hand side G
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nelem, &redsys.G));

    Mat BgT;
    PetscCall(MatCreateTranspose(redsys.Bg, &BgT));
    PetscCall(MatMult(BgT, redsys.g, redsys.G));
    PetscCall(VecScale(redsys.G, -1));

    return 1;
}

int DarcyStokes::CreateCoupledSystem(){

    int M1, N1, M2, N2;
    PetscCall(VecGetSize(redsys1->F, &M1)); 
    PetscCall(VecGetSize(redsys1->G, &N1));
    PetscCall(VecGetSize(redsys2->F, &M2)); 
    PetscCall(VecGetSize(redsys2->G, &N2));



    return 1;
}
