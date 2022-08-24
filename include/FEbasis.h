#ifndef FEBASIS_H_
#define FEBASIS_H_


#include <vector>
#include <set>
#include <array>
#include <valarray>
#include <algorithm>
#include <numeric>
#include <memory>

#include "lapacke.h"
#include "integral.h"
#include <assert.h>

using namespace std;

using vertex         = valarray<double>;
using vertexIndex    = valarray<int>;
using vertexSet      = set<vertex>;
using vertexIndexSet = set<vertexIndex>;

using namespace std;

typedef struct {


}

class MFE_Element{
    public:
        MFE_Element();

    private:
        // corners forming the element
        vertexSet corners_;

        vertexSet edgeMids_;
        vertexSet edgeUnitNormal_;
        vertexSet edgeUnitTang_;

        vertexSet diagNormal_; 



};

#endif
