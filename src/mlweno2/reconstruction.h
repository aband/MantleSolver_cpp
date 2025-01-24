#ifndef RECONSTMLWENO_H_
#define RECONSTMLWENO_H_

#include "stencilpolynomial.h"

class reconstruction{

    public:
        reconstruction() {};
        ~reconstruction() {};

        int prepare(const vector<int>& insize, const MeshInfo& mi);


    private:

        // Reconstruction stencil size
        vector<int> size {-1,-1};

        vector<stencilpolynomial> stencilPoly;

        // Map global index to local index
        map<int, int> globalTolocal;

        // Initialize 

};

#endif
