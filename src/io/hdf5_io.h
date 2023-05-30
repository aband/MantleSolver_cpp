#ifndef HDF5_IO_H_
#define HDF5_IO_H_

#include "hdf5.h"
#include "petsc.h"
#include <assert.h>

#include "util.h"

PetscErrorCode hdf5output(DM dmu, Vec* globalu);

PetscErrorCode phdf5Write(const MeshInfo& mi, DM dm, Vec * fullmesh);

#endif
