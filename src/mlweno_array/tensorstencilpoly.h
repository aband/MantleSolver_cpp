#ifndef TENSORSTENCILPOLY_H_
#define TENSORSTENCILPOLY_H_

#include "util.h"
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
        int setCoef(const MeshInfo& mi,
                    const int& gstartx, const int& gstarty);

    private:
        double * coef = nullptr;

        int order = 0;
        int sizex = 0;
        int sizey = 0;

        vertex center = {0.0,0.0};
        double h = 0.0;
};

#endif
