#ifndef TRANS_PARAM_H_
#define TRANS_PARAM_H_

#include "reconstMLWENO.h"

enum Location {LeftBndry, RightBndry, TopBndry, BottomBndry, Interior};

const std::array<double,2> advFunc(double u);

const std::array<double,2> dAdvFunc(double u);

const double diffFunc(double u);

const double dDiffFunc(double u);

double Distribution(const vertex& point,
                    const vector<double>& param);

bool left_boundary(const indice& globalCell,
                   const MeshInfo& mi);

bool right_boundary(const indice& globalCell,
                    const MeshInfo& mi);

bool top_boundary(const indice& globalCell,
                  const MeshInfo& mi);

bool bottom_boundary(const indice& globalCell,
                     const MeshInfo& mi);

bool interior(const indice& globalCell,
              const MeshInfo& mi);

// Assign location to mlwenouse object 
Location assignLocation(const indice& globalCell);

#endif
