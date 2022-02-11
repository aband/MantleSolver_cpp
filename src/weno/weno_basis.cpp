#include "../../include/weno.h"

double constfunc(valarray<double>&){
    return 1.0;
}

/*
 *The following member functions are invariant to different index, and
 *reconstruction order.
 */

void WenoBasisCoeff::SetStencilInfo(vector< valarray<int> >& si, 
                    vector<int>& o, vector<int>& r){
    order = o;
    stencilindex = si;
    range = r;
}

void WenoBasisCoeff::GetStencil(){

    stencil = new vector< vector< vector< valarray<double> > > >();
    stencilcenter = new vector< valarray<double> >();
    stencilh = new vector<double>();

    int totali = M+2*ghost;

    for (int j=range[2]; j<range[3]; j++){
    for (int i=range[0]; i<range[1]; i++){

        // Locate reference center point for each stencil
        vector< vector< valarray<double> > > stencilcell;
        valarray<double> scenter {0.0,0.0};

        for (auto & s: stencilindex){
            vector< valarray<double> > cellcorner;
            valarray<double> cellcenter {0.0,0.0};

            for (auto & c: cornerindex){
                valarray<double> corner = lmesh[(j+c[1]+s[1])*totali+(i+c[0]+s[0])];
                cellcenter += corner/(4.0); 
                cellcorner.push_back(corner);
            }

            scenter +=cellcenter/(double)(order[0]*order[1]);

            stencilcell.push_back(cellcorner);
        }

        double h = pow(NumIntegralFace(stencilcell.at(0),{0.0,0.0},1.0,constfunc),0.5);

        stencil->push_back(stencilcell);
        stencilcenter->push_back(scenter);  
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
