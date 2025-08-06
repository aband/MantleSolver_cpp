#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

#include "tensorstencilpoly.h"

// seriel version
class reconstruction {

    public:

        reconstruction()  {};
        ~reconstruction() {};

        int init(int sizex_sm, int sizey_sm,
                 int sizex_lg, int sizey_lg,
					  int order_sm, int order_lg,
					  const vector<indice>& sten_lg_pre,
                 const vector<indice>& sten_sm_pre,
                 const MeshInfo& mi, indice start);

        int setWgts(const vector<double>& stensigma_lg,
                    const vector<double>& stensigma_sm,
					  	  double h0);

        int extractsigma(const vector<double>& sigma_lg,
                         const vector<double>& sigam_sm);

        // Stencils
        vector<indice> sten_lg;
        vector<indice> sten_sm;
        int use_sten_const;

        // Linear weights
        vector<double> linwgts_lg;
        vector<double> linwgts_sm;
        double linwgts_const;

        // Nonlinear weights
        vector<double> nonlinwgts_lg;
        vector<double> nonlinwgts_sm;
        double nonlinwgts_const;

        int r_sm;
        int r_lg;
        int r_const = 1;

    private:

        vector<int> stencilnum;
        vector<double> nonlinwgts;

        double epsilon = 1e-4;
        int    s       = 1;
};

const int geteta(int r);

const bool validsten(int r, const MeshInfo& mi);

#endif
