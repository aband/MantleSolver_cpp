#ifndef LOCMATRIX_H_
#define LOCMATRIC_H_

#include <petsc.h>
#include "Hdivmixed.h"
#include "brmixed.h"
#include "util.h"
#include "shape.h"
#include "myFunc.h"

typedef struct{

  std::vector<double> A;
  std::vector<double> B;
  double C;
  std::vector<double> f;

} LocMat;

int CellAvePorosity(const MeshInfo& mi, 
                    PhysProperty * pp,
                    basis& basis_,
                    const valarray<double>& gwf,
                    const vector<vertex>& gpf);

int CellAvePorosity(const MeshInfo& mi, 
                    Phase * phase,
                    basis& basis_,
                    const valarray<double>& gwf,
                    const vector<vertex>& gpf);

// Stokes
int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 basis& basis_,
                 LocMat * loc,
                 PhysProperty * pp,
                 const valarray<double>& gwe, 
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 basis& basis_,
                 LocMat * loc,
                 Phase * phase,
                 const valarray<double>& gwe, 
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

// Darcy
int AssignLocMat(const MeshInfo& mi,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 LocMat * loc,
                 PhysProperty * pp,
                 const valarray<double>& gwe, 
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

int AssignLocMat(const MeshInfo& mi,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 LocMat * loc,
                 Phase * phase,
                 const valarray<double>& gwe, 
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

// Coupling term
int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 PhysProperty * pp,
                 double * k,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 Phase * phase,
                 double * k,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

// ===============================================================================

#ifdef COUPLED
#include "coupled.h"
int CellAvePorosity(const MeshInfo& mi, 
                    basis& basis_,
                    Phase * phase,
                    const indice& globalCell,
                    const valarray<double>& gwf,
                    const vector<vertex>& gpf,
                    MLWENO::MLWENOUse * mluse);

int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 basis& basis_,
                 LocMat * loc,
                 Phase * phase,
                 const indice& globalCell,
                 MLWENO::MLWENOUse * mluse,
                 const valarray<double>& gwe,
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

int AssignLocMat(const MeshInfo& mi,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 LocMat * loc,
                 Phase * phase,
                 const indice& globalCell,
                 MLWENO::MLWENOUse * mluse,
                 const valarray<double>& gwe,
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 double * k,
                 Phase * phase,
                 const indice& globalCell,
                 MLWENO::MLWENOUse * mluse,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf);

#endif



#endif
