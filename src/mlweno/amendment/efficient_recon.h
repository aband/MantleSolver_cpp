#ifndef EFFICIENT_RECON_H_
#define EFFICIENT_RECON_H_

#include "stencilpolynomial.h"
#include <unordered_map>

// Rearrange previously defined functions to make a more efficient reconstruction

class singlelevel{

    public:

        singlelevel() {};
        ~singlelevel() {};

        int prepare(int order, const vector<int>& insize, const MeshInfo& mi);

        int order = 0;

        int getsol(Tensor<double>& stencilsol, double ** localsol, 
                   const indice& stencilindex, const std::string& name) const;

        int updatesigma();

        // Point evaluation of each stencilpoly
        double evel(const vector<int>& stencilindex,
                    const vector<int>& baseindex,
                    const vertex& point) const;

        double eval(const vector<int>& index, 
                    const Tensor<double>& sol,
                    const vertex& point) const;

        // Print coefficients
        int printcoef(){for (int i=0;i<stencilPoly.getSize(); i++){stencilPoly(i).printCoef();} return 1;}

        int printCoef(int s){return stencilPoly(s).printCoef();}

    private:

        int getcenter(const vector<vector<vertex>>& cornerSet,
                      vector<vertex>& refcell,
                      vertex& center, double& h, double& area);
    
        Tensor<stencilpolynomial> stencilPoly;
        Tensor<double> sigma;

        int left, right, top, bottom;
};

#endif
