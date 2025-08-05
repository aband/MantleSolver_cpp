#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

#include "tensorstencilpoly.h"

// seriel version
class reconstruction {

    public:

        reconstruction()  {};
        ~reconstruction() {};

        int init();

        int setWgts(const vector<double>& stensigma, double h0);

        int extractsigma(const vector<double>& sigma_lg,
                         const vector<double>& sigam_sm);

        vector<indice> sten_lg;
        vector<indice> sten_sm;
        int use_sten_const;

        vector<double> linwgts;

        int count;

    private:

        vector<int> stencilnum;

        vector<double> nonlinwgts;

        vector<int> stenorder;

        double epsilon = 1e-4;
        int    s       = 1;

};

const int geteta(int r);

const bool validsten(int r, const MeshInfo& mi);

#endif
