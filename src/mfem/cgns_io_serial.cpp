#include "cgns_io_serial.h"

PetscErrorCode CgnsOutSerial(const MeshInfo& mi, const std::vector<double>& fullSol, 
                             int M, int N, basis& basis_, BRMixed& br, Hdivmixed& hdiv, 
                             PhysProperty * pp, int flag){

    double **x = new double*[N+1];
    double **y = new double*[N+1];

    double **vx = new double*[N];
    double **vy = new double*[N];

    for (int i=0; i<N+1; i++){
        x[i] = new double[M+1];
        y[i] = new double[M+1];
    }

    for (int j=0; j<N+1; j++){
        for (int i=0; i<M+1; i++){
            vertex vert = mi.lmesh[FlatIndic(M+1,i,j)];
            x[i][j] = vert[0];
            y[i][j] = vert[1];
        }
    }


    return PETSC_SUCCESS;
}
