#ifndef BOUNDARY_H_
#define BOUNDART_H_

#include <unordered_map>
#include <vector>
#include <array>
#include <valarray>
#include <algorithm>
#include <numeric>
#include <memory>

#include "lapacke.h"
#include "integral.h"
#include <assert.h>

using namespace std;

using point        = valarray<double>;
using point_index  = valarray<int>;
using index_set    = vector<point_index>;
using points_set   = vector<point>;
using stencil      = vector<points_set>;




#endif
