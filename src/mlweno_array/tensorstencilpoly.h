#ifndef TENSORSTENCILPOLY_H_
#define TENSORSTENCILPOLY_H_

#include "util.h"
#include "tensor.h"
#include <map>

class tensorstencilpoly {

    public:

        /*!
         * Highest polynomial order will be degree - 1
         */
        tensorstencilpoly() {};

        tensorstencilpoly(const int& order);

        tensorstencilpoly(const int& sizex, 
                          const int& sizey,
                          const int& order);

        ~tensorstencilpoly() {};

        // Compute stencil base polynomial coefficients
        int setCoef(const vector<vector<vertex>>& cornerSet, 
                    const vertex& center, const double& h);

        int setCoef(const vector<vector<vertex>>& cornerSet,
                    const double& x0, const double& y0, const double& h)
        {setCoef(cornerSet, {x0,y0}, h); return 1;};

    private:
        double * coef = nullptr;

        int order = 0;
        int sizex = 0;
        int sizey = 0;
}

#endif
