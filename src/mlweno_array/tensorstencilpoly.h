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

        ~tensorstencilpoly();

        // Compute stencil base polynomial coefficients
        int setCoef(const MeshInfo& mi, const int& gstartx, const int& gstarty);

        // Evaluate single stencil polynomials 
        double eval(const double& x, const double& y, const int& ncell) const;
        double eval(const double& x,  const double& y, 
                    const double& x0, const double& y0, 
                    const double& h,  const int& ncell) const;  

        int der(const int& derx, const int& dery,
                const double& x, const double& y,
                double* val);

        // Print stencil polynomial coefficients
        int printCoef();
        int printCoef(double* c, int n);

    private:
        double * coef = nullptr;

        int order = 0;
        int sizex = 0;
        int sizey = 0;

        vertex center  = {0.0,0.0};
        vector<vertex> refcell;
        double refarea = 0.0;
        double h = 0.0;
};

#endif
