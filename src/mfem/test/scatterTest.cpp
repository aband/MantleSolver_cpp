#include <iostream>
#include <petsc.h>

int main(int argc, char ** argv){

    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Vec Scatter Test \n"));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The code is running on %d processor(s) \n",size));

    // Create MPI Vector
	 Vec global; 
    int ownership = 2*(rank+1);
    VecCreateMPI(PETSC_COMM_WORLD, ownership, PETSC_DETERMINE, &global);
    VecSetUp(global);
    int low, high;
    VecGetOwnershipRange(global, &low, &high);

    for (int i=low; i<high; i++){
        double val = rank;
        VecSetValues(global, 1,&i, &val,INSERT_VALUES); 
    }

    VecAssemblyBegin(global);
    VecAssemblyEnd(global);

    int globalSize = 0;
    VecGetSize(global, &globalSize);

    VecView(global, PETSC_VIEWER_STDOUT_WORLD);

    // Create Scatter
    VecScatter scatter; 

    IS from, to;
    PetscInt *id_to  ;
    PetscInt *id_from;

    int nsize = globalSize;
    PetscMalloc1(nsize, &id_to);
    PetscMalloc1(nsize, &id_from);


    for (int i=0; i<nsize; i++){
        id_from[i] = i;
    }

    Vec destination;

    VecCreateSeq(PETSC_COMM_SELF, nsize, &destination);

    PetscScalar *values;

    ISCreateGeneral(PETSC_COMM_SELF, nsize, id_from, PETSC_COPY_VALUES, &from);
    ISCreateGeneral(PETSC_COMM_SELF, nsize, id_to, PETSC_COPY_VALUES, &to);

//    VecScatterCreate(global, from, destination, to, &scatter);
//    VecScatterBegin(scatter,global, destination, INSERT_VALUES, SCATTER_FORWARD);
//    VecScatterEnd(scatter,global, destination, INSERT_VALUES, SCATTER_FORWARD);

    VecScatterCreate(global, from, destination, NULL, &scatter);
    VecScatterBegin(scatter,global, destination, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter,global, destination, INSERT_VALUES, SCATTER_FORWARD);

	 VecGetArray(destination, &values);

    PetscPrintf(PETSC_COMM_SELF, " Current rank is %d, The scatterred value is %f \n", rank, values[0]);

    VecRestoreArray(destination, &values);

    VecView(destination, PETSC_VIEWER_STDOUT_SELF);

    VecDestroy(&destination);

    ISDestroy(&from);
    ISDestroy(&to);

    VecScatterDestroy(&scatter);

    PetscFinalize();
    return 0;
}
