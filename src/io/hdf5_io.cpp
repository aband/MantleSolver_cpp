// IO with hdf5 file format
#include "hdf5_io.h"

#define H5FILE_NAME "data.h5"
#define FAIL -1

PetscErrorCode hdf5output(DM dmu, Vec * globalu){

    PetscErrorCode   ierr;
    hid_t   file, dataset;       // file and dataset handles
    hid_t   datatype, dataspace; // handles
    hsize_t dimsf[2];            // dataset dimensions
    herr_t  status;

    Vec      fu, lu;
    PetscInt xs,ys,xm,ym,M,N,stencilwidth;
    double   **localu;

    PetscFunctionBeginUser;

    fu = *globalu;    

    ierr = DMDAGetCorners(dmu, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmu, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, &stencilwidth, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    int offsetx = M/4;
    int offsety = N/4;

    M = M-2*offsetx;
    N = N-2*offsety;

    double   *data = (double *)malloc(M*N*sizeof(double));

    // Get data
    ierr = DMGetLocalVector(dmu, &lu); CHKERRQ(ierr); 
    ierr = DMGlobalToLocalBegin(dmu,fu,INSERT_VALUES,lu); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmu,fu,INSERT_VALUES,lu); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmu,lu,&localu); CHKERRQ(ierr);

    for (int j=offsety; j<N+offsety; j++){
    for (int i=offsetx; i<M+offsetx; i++){
        data[(j-offsety)*M+(i-offsetx)] = localu[j][i];
    }}

    ierr = DMDAVecRestoreArrayRead(dmu,lu,&localu); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmu,&lu); CHKERRQ(ierr);

    /*
     *Create a new file using H5F_ACC_TRUNC access
     *H5F_ACC_TRUNC will cause the operation to overwrite
     *the previous file.
     */
    file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     *Describe the size of the array and create the data space for
     *fixed size dataset.
     */
    dimsf[0] = M; 
    dimsf[1] = N;
    dataspace = H5Screate_simple(2, dimsf, NULL);

    /*
     *Define datatype for the data in the file.
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE); 
    status   = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     *Create a new dataset within the file using defined
     *dataspace and datatype and default dataset creation properties.
     */
    dataset = H5Dcreate2(file, "Pressure", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     *Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    /*
     *Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);

    PetscFunctionReturn(0);
}

PetscErrorCode phdf5Write(const MeshInfo& mi, DM dm, Vec * fullmesh){

    /**
     * Be careful with hdf5 data format.
     * The first dimension of HDF5 is y direction.
     * The second dimension of HDF5 is x direction.
     */

    PetscFunctionBeginUser;

    hid_t    fid1;                                             /* HDF5 file IDs */
    hid_t    acc_tpl1;                                         /* File access templates */
    hid_t    xfer_plist;                                       /* Dataset transfer properties list */
    hid_t    sid1;                                             /* Dataspace ID */
    hid_t    file_dataspace;                                   /* File dataspace ID */
    hid_t    mem_dataspace;                                    /* memory dataspace ID */
    hid_t    dataset1, dataset2;                               /* Dataset ID */
    hsize_t  dims1[2] = {(hsize_t)mi.MPIglobalCellSize[1], 
                         (hsize_t)mi.MPIglobalCellSize[0]};    /* dataspace dim sizes */

    hsize_t start[2];            /* for hyperslab setting */
    hsize_t count[2], stride[2]; /* for hyperslab setting */

    herr_t ret; /* Generic return value */

    MPI_Comm comm = PETSC_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;



//! Set up hyperslab first

    PetscFunctionReturn(0);
}
