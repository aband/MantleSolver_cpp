#ifndef EFFICIENT_RECON_H_
#define EFFICIENT_RECON_H_

#include "stencilpolynomial.h"
#include <unordered_map>

class reconstruction{

    public:
				reconstruction() {};
				~reconstruction() {};


    private:
            vector<int> size {-1,-1};

            Tensor<stencilpolynomial> stencilPoly;

            int getcenter(const vector<vector<vertex>>& cornerSet,
                          vector<vertex>& refcell,
                          vertex& center, double& h, double& area);
}
