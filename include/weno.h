#ifndef WENO_BASIS_H_
#define WENO_BASIS_H_

#include <vector>
#include <valarray>
#include <algorithm>

#include "lapacke.h"
#include "integral.h"

/*
 *Define data types using in this code.
 */
using point        = valarray<double>;
using point_index  = valarray<int>;
using index_set    = vector<point_index>;
using cell_corners = vector<point>;
using stencil      = vector<cell_corners>;

using namespace std;

class WenoMesh{
    public:
        /*
         *Initialize class wenomesh with appropriate input.
         */
        WenoMesh(int M, int N, int g, vector< point >& lm) :
                 M(M), N(N), ghost(g), lmesh(lm) {};

    protected:
        int M, N;    // Local mesh size
        int ghost;   // Ghost layer thickness
        vector< point >& lmesh;
        /*
         *Loop inside a single element.
         */
        vector< point_index > corner_index {{0,0},{1,0},{1,1},{0,1}};
};

class WenoStencil : public WenoMesh{
    public:
        WenoStencil(int M, int N, int g, vector< point >& lm) : 
                    WenoMesh(M, N, g, lm) {};

        void GetStencilInfo(index_set& input_index_set, point_index& input_target_cell,
                            vector<int>& input_order);   

        void SetUpStencil();

        void PrintSingleStencil();

    protected:
        vector<int> polynomial_order;       // WENO restruction polynomial order, (xpow, ypow)
        /*
         *Index set used in calculation of the reconstruction is orderred in the
         *following. The members of this vector will be added to the starting index.
         */
        index_set index_set_stencil;     // Index set used in reconstruction
        /*
         *The cell targeted for reconstruction.
         */
        point_index target_cell;
        /*
         *Reference center point of the targeted cell.
         */
        point center;
        /*
         *Scale factor for each stencil.
         */
        double h;
        /*
         *Four corners of the target cell.
         */
        cell_corners target_cell_corners;
};

class WenoBasisCoeffStencil : public WenoStencil{
    public:
        WenoBasisCoeffStencil(int M, int N, int g, vector< point >& lm) : 
                              WenoStencil(M, N, g, lm) {};

        void CreateBasisCoeff();

        void PrintBasisCoeff();

    protected:
        double * wenobasiscoeff;
};

#endif
