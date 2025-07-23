#ifndef EFFICIENT_STENPOLY_H_
#define EFFICIENT_STENPOLY_H_

#include "polynomial.h"

class stenpoly{

    public:
        // Constructors
        stenpoly() {};

        stenpoly(const int& sizex,
                           const int& sizey); 

        stenpoly(const int& r) {stenpoly(r,r);};

        // Destructor
        ~stenpoly() {};

        vertex center;
        double h;

        int setCoef(const vector<vector<vertex>>& cornerSet,
                    const vertex& center, const double& scale);

    private:

        vector<int> size {-1,-1};

        Tensor<polynomial> tensorpoly;
};

#endif
