#include "../../include/weno.h"

double constfunc(valarray<double>&){
    return 1.0;
}

/*
 *The following member functions are invariant to different index, and
 *reconstruction order.
 */

void WenoBasisCoeff::GetStencil(){

    stencil = new vector< vector< vector< valarray<double> > > >();
    stencilcenter = new vector< valarray<double> >();
    stencilh = new vector<double>();

    for (int j=jbegin; j<jend; j++){
    for (int i=ibegin; i<iend; i++){

        // Locate reference center point for each stencil
        valarray<double> center {0.0,0.0};
        vector< vector< valarray<double> > > stencilcell;

        for (auto & s: stencilindex){
            vector< valarray<double> > cellcorner;

            for (auto & c: cornerindex){
                valarray<double> corner = lmesh[(j+c[1])*(M+2)+(i+c[0])];
                center += corner/(4.0*(double)order.size()); 
                cellcorner.push_back(corner);
            }

            stencilcell.push_back(cellcorner);
        }

        double h = pow(NumIntegralFace(stencilcell.at(0),center,1.0,constfunc),0.5);

        stencil->push_back(stencilcell);
        stencilcenter->push_back(center);  
        stencilh->push_back(h);
    }}

}

vector< vector<double> >& WenoBasisCoeff::CreateStencilIntegral(){

    vector< vector<double> > csi;

    return csi;
}

vector< vector<double> >& WenoBasisCoeff::CreateWenoBasisCoeff(){
    
    vector< vector<double> > wbc;

    return wbc;
}
