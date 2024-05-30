#include "preconst.h"

int SolScatAll(Vec * sol, Vec * g,
               Vec * destSol, Vec * destg){

    int totalDOF = 0;
    VecGetSize(*sol, &totalDOF);

    int bndryDOF = 0;
    VecGetSize(*g, &bndryDOF);

    VecScatter scatter_sol, scatter_g;

    IS from_sol, from_g, to_sol, to_g;
    PetscInt *id_from_sol;
    PetscInt *id_from_g;

    PetscMalloc1(totalDOF, &id_from_sol);
    PetscMalloc1(bndryDOF, &id_from_g);

    for (int i=0; i<totalDOF; i++){
        id_from_sol[i] = i;
    }

    for (int i=0; i<bndryDOF; i++){
        id_from_g[i] = i; 
    }

    VecCreateSeq(PETSC_COMM_SELF, totalDOF, destSol);
    VecCreateSeq(PETSC_COMM_SELF, bndryDOF, destg);

    ISCreateGeneral(PETSC_COMM_SELF, totalDOF, id_from_sol, 
                    PETSC_COPY_VALUES, &from_sol);
    ISCreateGeneral(PETSC_COMM_SELF, bndryDOF, id_from_g, 
                    PETSC_COPY_VALUES, &from_g);

    ISCreateGeneral(PETSC_COMM_SELF, totalDOF, id_from_sol, 
                    PETSC_COPY_VALUES, &to_sol);
    ISCreateGeneral(PETSC_COMM_SELF, bndryDOF, id_from_g, 
                    PETSC_COPY_VALUES, &to_g);

    VecScatterCreate(*sol, from_sol, *destSol, to_sol, &scatter_sol);
    VecScatterBegin(scatter_sol, *sol, *destSol, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter_sol, *sol, *destSol, INSERT_VALUES, SCATTER_FORWARD);

    VecScatterCreate(*g, from_g, *destg, to_g, &scatter_g);
    VecScatterBegin(scatter_g, *g, *destg, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter_g, *g, *destg, INSERT_VALUES, SCATTER_FORWARD);

    ISDestroy(&from_sol);
    ISDestroy(&from_g);

    VecScatterDestroy(&scatter_sol);
    VecScatterDestroy(&scatter_g);

    return 0;
}
