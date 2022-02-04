#ifndef WENO_BASIS_H_
#define WENO_BASIS_H_

#include <vector>
#include <valarray>
/*
 *#include "cblas.h"
 *#include "lapacke.h"
 */

using namespace std;

class wenobasis{
    public:
        // Constructors and class members
        int order;              // WENO restruction order
        valarray<double>& center; // Restruction reference center point
         /*
         *Index set used in calculation of the reconstruction is orderred in the
         *following. The members of this vector will be added to the starting index.
         */
        vector< valarray<int> >& index;     // Index set used in reconstruction

        wenobasis(valarray<double>& point, vector< valarray<int> >& o, int r) : 
                  center(point), index(o), order(r) {};

    private:
//        vector<double>& StencilIntegral;
};

#endif
