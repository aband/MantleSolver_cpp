#include "reconstruction.h"

const int geteta(int r){

    if (r==1){
        return 1;
    } else if (r==2){
        return 3;
    } else {
        return 4;
    }
}

// Return whether it is a valid stencil index
// In serial manner
const bool validsten(int sizex, int sizey, 
                     const MeshInfo& mi, indice stencil){

    int M = mi.MPIglobalCellSize[0] - sizex+1;
    int N = mi.MPIglobalCellSize[1] - sizey+1;

    if (stencil[0] < 0 || stencil[1] < 0 || stencil[0] > M-1 || stencil[1] > N-1 ){
        // it is not a valid stencil
        return false;
    } else {
        // it is a valid stencil
        return true;
    }
}

// =========================================================================

int reconstruction::init(int sizex_sm, int sizey_sm,
                         int sizex_lg, int sizey_lg,
                         int order_sm, int order_lg,
								 const vector<indice>& sten_lg_pre,
                         const vector<indice>& sten_sm_pre,
                         const MeshInfo& mi, indice start){

    sten_lg.clear();
    sten_sm.clear();

    allsize_lgx = mi.MPIglobalCellSize[0] - sizex_lg+1;
    allsize_lgy = mi.MPIglobalCellSize[1] - sizey_lg+1;

    allsize_smx = mi.MPIglobalCellSize[0] - sizex_sm+1;
    allsize_smy = mi.MPIglobalCellSize[1] - sizey_sm+1;

    for (auto it: sten_lg_pre){
				indice now = start + it;
        if (validsten(sizex_lg, sizey_lg, mi, now)){
				sten_lg.push_back(it);
				linwgts_lg.push_back(1);
				flat_sten_lg.push_back(now[1]*allsize_lgx + now[0]);
        }
    }

    for (auto it: sten_sm_pre){
				indice now = start + it;
        if (validsten(sizex_sm, sizey_sm, mi, now)){
				sten_sm.push_back(it);
            linwgts_sm.push_back(1);
				flat_sten_sm.push_back(now[1]*allsize_smx + now[0]);
        }
    }

    if (use_sten_const){
        linwgts_const = 0.0001;
    }

    // ===========================================================
    r_lg = order_lg + 1;
    r_sm = order_sm + 1;

    stensigma_lg.resize(linwgts_lg.size());
    stensigma_sm.resize(linwgts_sm.size());

    nonlinwgts_lg.resize(linwgts_lg.size());
    nonlinwgts_sm.resize(linwgts_sm.size());

    gstart = start;

    return 1;
}

int reconstruction::extractsigma(const vector<double>& sigma_lg, 
                                 const vector<double>& sigma_sm){

    //stensigma_lg.clear();
    //stensigma_sm.clear();

    for (int s=0; s<sten_lg.size(); s++){

//        indice stenid = sten_lg.at(s) + gstart;
//        stensigma_lg.at(s) = sigma_lg.at(stenid[1] * allsize_lgx + stenid[0]);
        stensigma_lg.at(s) = sigma_lg.at(flat_sten_lg.at(s));

    }

    for (int s=0; s<sten_sm.size(); s++){

//        indice stenid = sten_sm.at(s) + gstart;
//        stensigma_sm.at(s) = sigma_sm.at(stenid[1] * allsize_smx + stenid[0]);
        stensigma_sm.at(s) = sigma_sm.at(flat_sten_sm.at(s));
    }

    return 1;
}

int reconstruction::setWgts(double area){

    double sum = 0.0;

    //nonlinwgts_lg.clear(); 
    //nonlinwgts_sm.clear(); 

    for (int l=0; l<nonlinwgts_lg.size(); l++){
        nonlinwgts_lg.at(l) = linwgts_lg.at(l) / pow(stensigma_lg.at(l) + epsilon*area, s*r_lg + geteta(r_lg));
        sum += nonlinwgts_lg.at(l);
    }

    for (int l=0; l<nonlinwgts_sm.size(); l++){
        nonlinwgts_sm.at(l) = linwgts_sm.at(l) / pow(stensigma_sm.at(l) + epsilon*area, s*r_sm + geteta(r_sm));
        sum += nonlinwgts_sm.at(l);
    }

    if (use_sten_const){
        nonlinwgts_const = linwgts_const / pow(0.0 + epsilon*area, s*r_const + geteta(r_const));
		  sum += nonlinwgts_const;
    }

    for (int l=0; l<nonlinwgts_lg.size(); l++){
        nonlinwgts_lg.at(l) /= sum;
    }

    for (int l=0; l<nonlinwgts_sm.size(); l++){
        nonlinwgts_sm.at(l) /= sum;
    }

    nonlinwgts_const /= sum;

    return 1;
}

int reconstruction::eval(double ** localvals, const vector<vertex>& p, 
                         const vector<tensorstencilpoly>& sten_lg,
                         const vector<tensorstencilpoly>& sten_sm){

    elem_val.clear();
    elem_val.resize(p.size());

    double sum= 0.0;
    for (int c=0; c<p.size(); c++){
        sum = 0.0;

        for (int s=0; s<nonlinwgts_lg.size(); s++){
            sum += nonlinwgts_lg.at(s) * sten_lg.at(flat_sten_lg.at(s)).eval(localvals,p.at(c));
        } 

        for (int s=0; s<nonlinwgts_sm.size(); s++){

            sum += nonlinwgts_sm.at(s) * sten_sm.at(flat_sten_sm.at(s)).eval(localvals,p.at(c));
       }

        if (use_sten_const){
            sum += nonlinwgts_const * localvals[gstart[1]][gstart[0]];
        }

        elem_val.at(c) = sum;
    }

    return 1;
}

// ===========================================================
int reconstruction::printinfo(){

    cout << "Number of " << sten_lg.size() <<  " large stencil of order : " << r_lg  << " is used." << endl;

    for (int s=0; s<sten_lg.size(); s++ ){
        cout << "At stencil : (" << sten_lg.at(s)[0] + gstart[0] << ", " << 
                                    sten_lg.at(s)[1] + gstart[1] << ") " << 
												nonlinwgts_lg.at(s) << "  " << endl;
    }cout << endl;

    cout << "Number of " << sten_sm.size() <<  " small stencil of order : " << r_sm  << " is used." << endl;

    for (int s=0; s<sten_sm.size(); s++ ){
        cout << "At stencil : (" << sten_sm.at(s)[0] + gstart[0] << ", " << 
                                    sten_sm.at(s)[1] + gstart[1] << ") : " << 
												nonlinwgts_sm.at(s) << "  " << endl;
    }cout << endl;

    if (use_sten_const){
        cout << "The constant stencil : " << nonlinwgts_const << endl;
      
    }cout << endl;

    return 1;
}

double reconstruction::efforder(){

    double sum = 0;

    for (int s=0; s<sten_lg.size(); s++ ){
        sum += nonlinwgts_lg.at(s) * r_lg;
    }
  
    for (int s=0; s<sten_sm.size(); s++ ){
        sum += nonlinwgts_sm.at(s) * r_sm;
    }

    if (use_sten_const){
        sum += nonlinwgts_const * r_const;
    }

    return sum;
}
