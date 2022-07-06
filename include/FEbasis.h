#ifndef FEBASIS_H_
#define FEBASIS_H_


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

using namespace std;

typedef struct {


}

class Element{
    public:
        Element();

    private:
        vector<point> corner;

        vector<point> mid;
        vector<point> unit_normal;
        vector<point> unit_tang;

        

};

#endif
