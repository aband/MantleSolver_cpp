#include "cgns_io_serial.h"

PetscErrorCode CgnsOutSerial(const MeshInfo& mi, const std::vector<double>& fullSol, 
                             int M, int N, basis& basis_, BRMixed& br, Hdivmixed& hdiv, 
                             PhysProperty * pp, int flag){

    int index_file, icelldim, iphysdim, index_base;
    int index_zone, index_coord;
    char basename[33],zonename[33];

    double **x = new double*[N+1];
    double **y = new double*[N+1];

    double **vx = new double*[N];
    double **vy = new double*[N];

    for (int i=0; i<N+1; i++){
        x[i] = new double[M+1];
        y[i] = new double[M+1];
    }

    for (int i=0; i<N; i++){
        vx[i] = new double[M];
        vy[i] = new double[M];
    }

    // Assign grid
    for (int j=0; j<N+1; j++){
        for (int i=0; i<M+1; i++){
            vertex vert = mi.lmesh[FlatIndic(M+1,i,j)];
            x[j][i] = vert[0];
            y[j][i] = vert[1];
        }
    }

    if (cg_open("grid.cgns", CG_MODE_WRITE, &index_file)) cg_error_exit();

    strcpy(basename,"Base");

    icelldim = 2;
    iphysdim = 2;

    cg_base_write(index_file, basename, icelldim, iphysdim, &index_base);


    return PETSC_SUCCESS;
}
