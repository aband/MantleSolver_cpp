#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <petscdmda.h>
#include "mesh.h"

PetscErrorCode PrintFullMesh(DM dmMesh, Vec * fullmesh);

/**
 * Specified to output pressure data.
 * The output file name is pre determined as pressure.data.
 */
PetscErrorCode DrawPressure(DM dmu, Vec * globalu);

/**
 * Output data in a unprocessed raw way.
 * Almost the same as function DrawPressure expect filename is not pre determined.
 */
PetscErrorCode PlainOutput(DM dmu, Vec * globalu, char* filename);

/**
 * Output coordiante data in an unprocessed raw way.
 * Two files will be created, gridX and gridY.
 */
PetscErrorCode PlainMeshOutput(DM dmMesh, Vec * fullmesh);

/**
 * Output (Jacobian) matrix element pattern.
 */
PetscErrorCode DrawMat(Mat V, const char * myfile);

PetscErrorCode WriteMat(Mat V, const char * myfile);

PetscErrorCode WriteVec(Vec g, const char * myfile);

/**
 * Output with native mpi I/O format
 */
PetscErrorCode MPIIO();

#endif
