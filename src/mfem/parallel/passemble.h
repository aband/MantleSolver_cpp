#ifndef PASSEMBLE_H_
#define PASSEMBLE_H_

#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "myFunc.h"
#include "assemble.h"
#include "locmat.h"
#include "bndry.h"
#include "solve.h"
#include "shape.h"

PetscErrorCode ParallelAssembleTest();

template <typename T>
inline int CreateRefMap(T& funcSp, int * refArray, 
                        const MeshInfo& mi, int * bndryDOFEssen){

    // Each processor has to create its own mapping
    // Control Essential dof only

    // !!!!!! Caution !!!!!!!
    // This function has not been finished
    // Cannot do natural boundary condition yet

    int bndryIndex = 0;
    int intrIndex  = 0;

    for (int dof=0; dof<funcSp.getDOF(); dof++){

        if (funcSp.onBndry(mi, dof)){

            refArray[dof] = intrIndex;
            intrIndex ++;
        }else {
            refArray[dof] = bndryIndex;
            bndryIndex ++;
        }

    }

    *bndryDOFEssen = intrIndex;

    return 0;  
};

// ! A completed version of assigning boundary essential conditions
template <typename T>
inline int CreateRefMap(T& funcSp, 
                        const MeshInfo& mi, 
                        bool (*EssenBndry)(const MeshInfo&, T&, const int&),
                        int * refArray,
                        unordered_map<int, int>& refMapNatur,
                        int * EssenDOFCount,
                        int * NaturDOFCount){

    int essenCount = 0; 
    int naturCount = 0;
    int interCount = 0;

    // ! loop through all dofs 
    // Natural boundary dofs are different from Essential boundary dofs
    // For the fact that natrual boundary dofs participate in the left hand side matrix
    // and also right hand side vector
    // Hence it requires two different index system.
    for (int dof=0; dof<funcSp.getDOF(); dof++){
        if (funcSp.onBndry(mi, dof)){

            if (EssenBndry(mi, funcSp, dof)){
                refArray[dof] = essenCount;
                essenCount ++;
            } else {
                refArray[dof] = interCount;
                interCount ++;

                refMapNatur.insert(std::make_pair<int, int>(dof, naturCount));
                naturCount ++;
            }

        } else {
                refArray[dof] = interCount;
                interCount ++;
        }

    }

    *EssenDOFCount = essenCount;
    *NaturDOFCount = naturCount;

    return 0;
};

PetscErrorCode ParallelMatrixAssemble(const MeshInfo& mi,
                                      basis& basis_,
                                      PhysProperty * pp,
                                      const bndryVal& bndryEssenStokes,
                                      ReducedSys * redsysStokes,
                                      const bndryVal& bndryEssenDarcy,
                                      ReducedSys * redsysDarcy,
                                      Mat * K,
                                      BRMixed& br_,
                                      Hdivmixed& hdiv_,
                                      int * refArrayStokes,
                                      int * refArrayDarcy,
                                      const int& bndryDOFStokes,
                                      const int& bndryDOFDarcy);

PetscErrorCode ParallelMatrixAssemble(const MeshInfo& mi,
                                      basis& basis_,
                                      Phase * phase,
                                      const bndryVal& bndryEssenStokes,
                                      ReducedSys * redsysStokes,
                                      const bndryVal& bndryEssenDarcy,
                                      ReducedSys * redsysDarcy,
                                      Mat * K,
                                      BRMixed& br_,
                                      Hdivmixed& hdiv_,
                                      int * refArrayStokes,
                                      int * refArrayDarcy,
                                      const int& bndryDOFStokes,
                                      const int& bndryDOFDarcy);

// ======= Inline functions =====================

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

inline int PrepareReducedSys(ReducedSys * redsys, 
                             int reducedDOF, int bndrySize, int totalElem,
                             int Adnz, int Aonz,
                             int Bdnz, int Bonz){

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
    int m, n; 
    PetscCall(MatGetOwnershipRange(redsys->Bg, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys->g));
    PetscCall(VecSetUp(redsys->g));

    PetscCall(MatGetOwnershipRange(redsys->B, &m, &n));
    PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n-m, PETSC_DETERMINE, &redsys->source));
    PetscCall(VecSetUp(redsys->source));

    return 0;
}

inline int AssembleReducedSys(ReducedSys * redsys){

    PetscCall(MatAssemblyBegin(redsys->M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(redsys->Kg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->Kg, MAT_FINAL_ASSEMBLY));

    PetscCall(MatAssemblyBegin(redsys->B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(redsys->Bg, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->Bg, MAT_FINAL_ASSEMBLY));

    PetscCall(MatAssemblyBegin(redsys->C, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(redsys->C, MAT_FINAL_ASSEMBLY));

    PetscCall(VecAssemblyBegin(redsys->g));
    PetscCall(VecAssemblyEnd(redsys->g));

    PetscCall(VecAssemblyBegin(redsys->source));
    PetscCall(VecAssemblyEnd(redsys->source));

    return 0;
}

// Used for elements on the boundary (has dofs on the boundary)
template <typename T> 
inline int AssignLocRedSys(ReducedSys * redsys,
                           LocMat * loc,
                           int * ref,
                           const MeshInfo& mi,
                           const bndryVal& bndryEssen,
                           const indice& global,
                           shape<T>& funcSp){

    // ! Get global index of local dofs 
    const std::vector<int> elemDofs = funcSp.LocalToGlobal(mi, global);

    const int idxn = FlatIndic(mi, global);

    for (int row=0; row<elemDofs.size(); row++){
        // View it as the row index
        const int idxm = ref[elemDofs.at(row)];
        const double valB = loc->B.at(row);

        if (funcSp.onBndry(mi, elemDofs.at(row))){
            // This dof is on the boundary
            // Should be assign to Bg
            PetscCall(MatSetValues(redsys->Bg, 1, &idxm, 1, &idxn, &valB, 
                                   ADD_VALUES));

            // At the same time insert essential boundary value to rhs vector
            auto itFind = bndryEssen.find(elemDofs.at(row));
            if (itFind != bndryEssen.end()){
                const bndryInfo& tmp = bndryEssen.at(elemDofs.at(row));
                PetscCall(VecSetValues(redsys->g, 1, &idxm, &tmp.val, INSERT_VALUES));
            }
        } else {
            // This dof is not on the boundary
            // Should be assigned to B instead
            PetscCall(MatSetValues(redsys->B , 1, &idxm, 1, &idxn, &valB, 
                                   ADD_VALUES));

            // This dof is not on the boundary
            // This dof will contribute to source term
            double vals = loc->f.at(row);
            PetscCall(VecSetValues(redsys->source, 1, &idxm, &vals, ADD_VALUES));

            for (int col=0; col<elemDofs.size(); col++){
                const int cidxn  = ref[elemDofs.at(col)];
                const double val = loc->A.at(row+col*elemDofs.size()); 

                if (funcSp.onBndry(mi,elemDofs.at(col))){
                    // It is a non bndry dof - bndry dof interaction
                    // Val assigned to M
                    PetscCall(MatSetValues(redsys->Kg, 1, &idxm, 1, &cidxn, 
                                           &val, ADD_VALUES));
                } else {
                    // It is a non bndry dof - non bndry dof interaction
                    // Val assigned to M
                    PetscCall(MatSetValues(redsys->M, 1, &idxm, 1, &cidxn, 
                                           &val, ADD_VALUES));
                }
            }
        }
    }

    return 0;
}

// Used for elements on the boundary (Separating essential and natural boundary conditions)
template <typename T>
inline int AssignLocRedSys(ReducedSys * redsys, 
                           LocMat * loc,
                           int * ref,
                           const unordered_map<int, int>& refMapNatur,
                           const MeshInfo* mi,
                           const bndryVal& bndryEssen,
                           const indice& global,
                           bool (*EssenBndry)(const MeshInfo& mi, T&, const int&)
                           T& funcSp){

    // ! Get global index of local dofs 
    const std::vector<int> elemDofs = funcSp.LocalToGlobal(mi, global);

    const int idxn = FlatIndic(mi, global);

    for (int row=0; row<elemDofs.size(); row++){
        // View it as the row index
        // There are three meanings of values in ref array
         
        const int idxm = ref[elemDofs.at(row)];
        const double valB = loc->B.at(row);
    
        if (funcSp.onBndry(mi, elemDofs.at(row))){
            // The dof is on the boundary
            // Need further indenfication whether it is essential or natural
            if (EssenBndry(mi, T, elemDofs.at(row))){
                // If it if essential boundary dof it goes to right hand side vector g
                PetscCall(MatSetValues(redsys->Bg, 1, &idxm, 1, &idxn, &valB, 
                                       ADD_VALUES));


            }


        } else {

        }

    }

}

// Used for interior elements (no need to identify boundary dofs)
// Set multiple values at the same time
// need const int id array
// Stokes and Darcy part are separated
inline int AssignLocRedSys(ReducedSys * redsys,
                           LocMat * loc,
                           int * ref,
                           const MeshInfo& mi,
                           const indice& global,
                           Hdivmixed& hdiv_){

    const int idxm = FlatIndic(mi, global);

    std::array<int, 8> locDof = hdiv_.LocalToGlobal(mi, global);

    const int IS[8] = {ref[locDof[0]], ref[locDof[1]], 
                       ref[locDof[2]], ref[locDof[3]],
                       ref[locDof[4]], ref[locDof[5]], 
                       ref[locDof[6]], ref[locDof[7]]};

    const double locB[8] = {loc->B.at(0), loc->B.at(1),
                            loc->B.at(2), loc->B.at(3),
                            loc->B.at(4), loc->B.at(5),
                            loc->B.at(6), loc->B.at(7)};

    // Only fill B and M matrix
    PetscCall(MatSetValues(redsys->B, 8, IS, 1, &idxm, locB, ADD_VALUES));

    for (unsigned int l=0; l<8; l++){
        const double locA[8] = {loc->A.at(0+l*8), loc->A.at(1+l*8),
                                loc->A.at(2+l*8), loc->A.at(3+l*8),
                                loc->A.at(4+l*8), loc->A.at(5+l*8),
                                loc->A.at(6+l*8), loc->A.at(7+l*8)};
        const int Aidxm = IS[l];
        PetscCall(MatSetValues(redsys->M,1,&Aidxm,8,IS,locA,ADD_VALUES));
    }

    // Fill Source vector
    const double locf[8] = {loc->f.at(0), loc->f.at(1),
                            loc->f.at(2), loc->f.at(3),
                            loc->f.at(4), loc->f.at(5),
                            loc->f.at(6), loc->f.at(7)};

    PetscCall(VecSetValues(redsys->source, 8, IS, locf, ADD_VALUES));

    return 0;
}

inline int AssignLocRedSys(ReducedSys * redsys,
                           LocMat * loc,
                           int * ref,
                           const MeshInfo& mi,
                           const indice& global,
                           BRMixed& br_){

    const int idxm = FlatIndic(mi, global);

    std::array<int, 12> locDof = br_.LocalToGlobal(mi, global);

    const int IS[12] = {ref[locDof[0]], ref[locDof[1]], 
                        ref[locDof[2]], ref[locDof[3]],
                        ref[locDof[4]], ref[locDof[5]], 
                        ref[locDof[6]], ref[locDof[7]],
                        ref[locDof[8]], ref[locDof[9]],
                        ref[locDof[10]], ref[locDof[11]]};

    const double locB[12] = {loc->B.at(0), loc->B.at(1),
                             loc->B.at(2), loc->B.at(3),
                             loc->B.at(4), loc->B.at(5),
                             loc->B.at(6), loc->B.at(7),
                             loc->B.at(8), loc->B.at(9),
                             loc->B.at(10), loc->B.at(11)};

    // Only fill B and M matrix
    PetscCall(MatSetValues(redsys->B, 12, IS, 1, &idxm, locB, ADD_VALUES));

    for (unsigned int l=0; l<12; l++){
        const double locA[12] = {loc->A.at(0+l*12), loc->A.at(1+l*12),
                                 loc->A.at(2+l*12), loc->A.at(3+l*12),
                                 loc->A.at(4+l*12), loc->A.at(5+l*12),
                                 loc->A.at(6+l*12), loc->A.at(7+l*12),
                                 loc->A.at(8+l*12), loc->A.at(9+l*12),
                                 loc->A.at(10+l*12), loc->A.at(11+l*12)};
        const int Aidxm = IS[l];
        PetscCall(MatSetValues(redsys->M,1,&Aidxm,12,IS,locA,ADD_VALUES));
    }

    // Only contribute to source vector not g vector

    const double locf[12] = {loc->f.at(0), loc->f.at(1),
                             loc->f.at(2), loc->f.at(3),
                             loc->f.at(4), loc->f.at(5),
                             loc->f.at(6), loc->f.at(7),
                             loc->f.at(8), loc->f.at(9),
                             loc->f.at(10), loc->f.at(11)};

    PetscCall(VecSetValues(redsys->source, 12, IS, locf, ADD_VALUES));

    return 0;
}

#endif
