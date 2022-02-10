#ifndef WENO_BASIS_H_
#define WENO_BASIS_H_

#include <vector>
#include <valarray>
#include "integral.h"
/*
 *#include "cblas.h"
 *#include "lapacke.h"
 */

using namespace std;

class WenoMesh{
    public:
        int M, N;
        int ghost;
        vector< valarray<double> >& lmesh;
        /*
         *Loop inside a single element.
         */
        vector< valarray<double> > index {{0.0,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}};
        /*
         *Initialize class wenomesh with appropriate input.
         */
        WenoMesh(int M, int N, vector <valarray<double> >& lm) :
                 M(M), N(N), lmesh(lm) {};
};

class WenoBasisCoeff{
    public:
        // Constructor
        WenoBasisCoeff(valarray<double>& point, vector< valarray<int> >& o, int r) : 
                       center(point), index(o), order(r) {};

        void CreateStencilIntegral();
        void CreateWenoBasisCoeff();

    private:
        int order;                // WENO restruction order
        valarray<double>& center; // Restruction reference center point
         /*
         *Index set used in calculation of the reconstruction is orderred in the
         *following. The members of this vector will be added to the starting index.
         */
        vector< valarray<int> >& index;     // Index set used in reconstruction

        vector<double> * stencilintegral;
        vector<double> * coeff;
};

/*
 *class SmoothnessIndicator{
 *
 *};
 *
 *class WenoReconstruction: public WenoMesh, public WenoBasisCoeff, public SmoothnessIndicator{
 *
 *};
 */


#endif
