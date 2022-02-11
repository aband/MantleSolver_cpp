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
        /*
         *Initialize class wenomesh with appropriate input.
         */
        WenoMesh(int M, int N, int g, vector <valarray<double> >& lm) :
                 M(M), N(N), ghost(g), lmesh(lm) {};

    protected:
        int M, N;    // Local mesh size
        int ghost;   // Ghost layer thickness
        vector< valarray<double> >& lmesh;
        /*
         *Loop inside a single element.
         */
        vector< valarray<int> > cornerindex {{0,0},{1,0},{1,1},{0,1}};
};

class WenoBasisCoeff : public WenoMesh{
    public:
        // Constructor
        WenoBasisCoeff(int M, int N, int g, vector <valarray<double> >& lm) :
                       WenoMesh(M, N, g, lm) {};

        void SetStencilInfo(vector< valarray<int> >& stencilindex, 
                            vector<int>& order, vector<int>& range);

        // Get stencil information. Need to be called first!
        void GetStencil();

        // Create Integral for a given stencil
        vector< vector<double> >& CreateStencilIntegral();

        // Create weno basis polynomial coefficients based on the given stencil
        vector< vector<double> >& CreateWenoBasisCoeff();



    private:
        vector<int> order;       // WENO restruction polynomial order, (xpow, ypow)
         /*
         *Index set used in calculation of the reconstruction is orderred in the
         *following. The members of this vector will be added to the starting index.
         */
        vector< valarray<int> > stencilindex;     // Index set used in reconstruction

        /*
         *range of index that include start index and end index of i and j.
         *ibegin,iend,jbegin,jend.
         */
        vector<int> range;

        /*
         *Stencil contains four corners of each cell within the defined
         *stencil. Center point is defined by averaging out all points
         *in stencil.
         */
        vector< vector< vector< valarray<double> > > > * stencil;
        vector< valarray<double> > * stencilcenter;
        vector<double> * stencilh;
};

/*
 *class WenoReconstruction : public WenoBasisCoeff{
 *    public:
 *        WenoReconstruction();
 *
 *
 *}
 */

#endif
