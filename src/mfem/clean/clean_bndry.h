#ifndef CLEAN_BNDRY_H_
#define CLEAN_BNDRY_H_

#include "myFunc.h"
#include "util.h"
#include "Hdivmixed.h"
#include "brmixed.h"

// The boundary value data structure contains
// 1. global index of degree of freedom and global index of element
// 2. a pair object pairing local degree of freedom and value
struct bndryInfo{
    int        localDOF;
    double     val;  
    indice     globalElem;
};

using bndryVal = std::unordered_map<int, bndryInfo>;

int MarkBndryDOFStokes(bndryVal& essenVal,
                       bndryVal& naturVal,
                       const MeshInfo& mi,
                       basis& basis_,
                       BRMixed& br_,
                       PhysProperty * pp,
                       const std::vector<double>& paramter);

int MarkBndryDOFDarcy(bndryVal& essenVal,
                      bndryVal& naturVal,
                      const MeshInfo& mi,
                      basis& basis_,
                      Hdivmixed& hdiv_,
                      PhysProperty * pp,
                      const std::vector<double>& paramter);

#endif
