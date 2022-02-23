#ifndef WENO_BASIS_H_
#define WENO_BASIS_H_

#include <vector>
#include <valarray>
#include <algorithm>
#include <numeric>

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

/*
 *The number of variable will not always be zero.
 */
using solution = double;

using namespace std;

class WenoMesh{
    public:
        /*
         *Initialize class wenomesh with appropriate input.
         */
        WenoMesh(int M, int N, int g, vector< point >& lm, solution**& lsol) :
                 M(M), N(N), ghost(g), lmesh(lm), lsol(lsol) {};

        int M, N;    // Local mesh size
        int ghost;   // Ghost layer thickness
        vector< point >& lmesh;
        solution**& lsol;
};

class WenoStencil{
    public:
        WenoStencil(index_set& input_index_set, point_index& input_target_cell,
                    vector<int>& input_order):
                    polynomial_order(input_order), 
                    index_set_stencil(input_index_set), 
                    target_cell(input_target_cell) {};

        void SetUpStencil(const WenoMesh*& wm);

        void PrintSingleStencil();

        /*
         *Reference center point of the targeted cell.
         */
        point center;
        /*
         *Scale factor for each stencil.
         */
        double h;
        /*
         *Index set used in calculation of the reconstruction is orderred in the
         *following. The members of this vector will be added to the starting index.
         */
        index_set index_set_stencil;     // Index set used in reconstruction

    protected:

        vector<int> polynomial_order;       // WENO restruction polynomial order, (xpow, ypow)
        /*
         *The cell targeted for reconstruction.
         */
        point_index target_cell;
        /*
         *Four corners of the target cell.
         */
        cell_corners target_cell_corners;

        /*
         *Loop inside a single element.
         */
        vector< point_index > corner_index {{0,0},{1,0},{1,1},{0,1}};
};

class WenoPrepare : public WenoStencil{
    public:
        WenoPrepare(index_set& input_index_set, point_index& input_target_cell,
                    vector<int>& input_order) : 
                    WenoStencil(input_index_set,input_target_cell,input_order) {};

        void CreateBasisCoeff(const WenoMesh*& wm);

        void CreateSmoothnessIndicator(const WenoMesh*& wm, double eta, double Theta);

        void PrintBasisCoeff();

        /*
         *Define parameters
         */

        double omega = 0.0;
        double * wenobasiscoeff;

        double sigma = 0.0;
    private:

        double epsilon_0 = 0.5;
};

// Define independent functions
solution WenoReconstStencil(vector<int>& order, point_index& target, point target_point,
                            WenoPrepare*& wp,const WenoMesh*& wm);

solution WenoPointReconst(index_set& StencilLarge, vector<index_set>& StencilSmall,const WenoMesh*& wm,
                          point_index& target, vector<int>& Sorder, vector<int>& Lorder,
                          point target_point);


#endif
