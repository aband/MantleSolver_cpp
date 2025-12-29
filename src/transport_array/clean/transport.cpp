#include "transport.h"

int TransportVariable::CreateReconstruction(const MeshInfo& mi, 
                                            int sizelgx, int sizelgy, int orderlg,
                                            int sizesmx, int sizesmy, int ordersm){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    int Mlg = M-sizelgx+1;
    int Nlg = N-sizelgy+1;

    int Msm = M-sizesmx+1;
    int Nsm = N-sizesmy+1;

    stenlg.resize(Mlg*Nlg);

    for (int j=0; j<Nlg; j++){
    for (int i=0; i<Mlg; i++){
        int s = j*Mlg+i;
        stenlg.at(s) = tensorstencilpoly(orderlg, sizelgx, sizelgy);
        stenlg.at(s).setCoef(mi,i,j);
		  stenlg.at(s).setSigma();
		  stenlg.at(s).startx = i;
		  stenlg.at(s).starty = j;
    }}

    stensm.resize(Msm*Nsm);

    for (int j=0; j<Nsm; j++){
    for (int i=0; i<Msm; i++){
        int s = j*Msm + i;
        stensm.at(s) = tensorstencilpoly(ordersm, sizesmx, sizesmy);
        stensm.at(s).setCoef(mi,i,j);
		  stensm.at(s).setSigma();
		  stensm.at(s).startx = i;
		  stensm.at(s).starty = j;
    }}


    return 1;
}
