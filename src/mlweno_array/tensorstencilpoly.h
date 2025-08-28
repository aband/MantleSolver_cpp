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

        tensorstencilpoly(const int& inorder,
                          const int& insizex, 
                          const int& insizey);

        ~tensorstencilpoly();

        // Compute stencil base polynomial coefficients
        int setCoef(const MeshInfo& mi, const int& gstartx, const int& gstarty);

        // Compute sigma tensor base
        int setSigma();

        int setSigma(const MeshInfo& mi, indice local);

        int setSigma_p(double ** localvals);

        // Evaluate single stencil polynomials 
        double eval(const double& x, const double& y, const int& ncell) const;
        double eval(const double& x,  const double& y, 
                    const double& x0, const double& y0, 
                    const double& h,  const int& ncell) const;  

        double eval(double ** localsol,
                    const int& startx, const int& starty,
                    const double& x, const double& y) const;

        double eval(double ** localsol, const vertex& p) const {return eval(localsol, startx, starty, p[0], p[1]);};

        int der(int derX, int derY, int ncell, double * dp, double x, double y);
        int der(int derX, int derY, int ncell, double * dp, double x, double y, double x0, double y0, double scale);

        double sigma(double ** localsol, const int& startx, const int& starty);

        double sigma(double ** localsol){return sigma(localsol, startx, starty);};

        int sigma(double ** localsol, vector<double>& locsigma);

        int startx;
        int starty;

        // Print stencil polynomial coefficients
        int printCoef();
        int printCoef(double* c, int n);

        int printcollapseCoef(double ** localsol) const;

        int printSigmaBase();

        int printStencilSol(double ** localsol) const;

    private:
        double * coef = nullptr;

        double * sigmabase = nullptr;

        unordered_map<int, double *> sigmabasetarget;

        int order = 0;
        int sizex = 0;
        int sizey = 0;

        vertex center  = {0.0,0.0};
        vector<vertex> refcell;
        double refarea = 0.0;
        double h = 0.0;

};

#endif
