#include "../../include/weno.h"

double constfunc(valarray<double>&){
    return 1.0;
}

void CreateWenoBasis(wenomesh &wm){

    vector<wenobasis> weno2;
    vector<wenobasis> weno3;
    vector<wenobasis> weno4H;
    vector<wenobasis> weno4V;

    int M = wm.M;
    int N = wm.N;
    int ghost = wm.ghost;

    /*
     *2nd order weno reconstruction
     */
    for (size_t j=ghost-1; j<N+ghost+1; j++){
        for (size_t i=ghost-1; i<M+ghost+1; i++){
            valarray<double> center = wm.lmesh[j*(M+2)+i]; 

            vector< valarray<double> > corner;
            for (auto & p: wm.index){
                int newi = i+p[0];
                int newj = j+p[1];
                corner.push_back(center + wm.lmesh[newj*(M+2)+newi]);
            }

            double h= pow(NumIntegralFace(corner, center, 1.0, constfunc),0.5); 
            
            vector< valarray<int> > order {{-1,-1},{0,-1},{0,0},{-1,0}};
            wenobasis wb = wenobasis(center, order, 2);
            //wb.CreateStencilIntegral();
            //wb.CreateWenoBasisCoeff();
            weno2.push_back(wb);
        }
    }


}
