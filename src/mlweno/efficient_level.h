#ifndef EFFICIENT_LEVEL_H_
#define EFFICIENT_LEVEL_H_

#include "stenpoly.h"
#include <unordered_map>

class reconlevel{

    public:
        reconlevel() {};
        ~reconlevel() {};

        int prepare(const vector<int>& insize, const MeshInfo& mi);

    private:

        vector<int> size {-1,-1};
        Tensor<stencilpoly> sp;

        int left, right, bottom, top;
}

#endif
