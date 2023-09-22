#include "driver.h"

double func(const vertex& point, const vector<double>& param){
	 if (point[0]<param[0]){
//		  return point[0]*point[0]+point[1]*point[1];
	     return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1);
//        return sin(point[0] + point[1] + 0.1);
	 } else {
//		  return point[0]*point[0]*point[1]*point[1] + 1.0;
	     return sin(point[0]*3.0+0.5)+cos(point[1]/2.0-0.2) + pow(point[0]+0.1,3)*(point[1]+1) + 10;
//        return sin(point[0] + point[1] + 0.1) + 10;
    }

    // Infinitly smooth test case
    //return sin(point[0] + point[1] + 0.1);

    // Quadratic test case
    //return point[0]*point[0] + point[1]*point[1];

    // Linear test case
    //return point[0] + point[1];

    // Constant test case
    //return 0.5;
}

int main(int argc, char **argv){

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    Driver * driverPtr = new Driver();

    driverPtr->Prepare();

    driverPtr->CellAveragedInit(func);

    PetscFinalize();

    return 0;
}
