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

#endif
