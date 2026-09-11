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

    // clear linear weights
    linwgts_lg.clear();
    linwgts_sm.clear();
    flat_sten_lg.clear();
    flat_sten_sm.clear();

    sten_lg.clear();
    sten_sm.clear();

    stensizelgx = sizex_lg;
    stensizelgy = sizey_lg;

    stensizesmx = sizex_sm;
    stensizesmy = sizey_sm;

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
        linwgts_const = 0.00001;
    } else {
        linwgts_const = 0.000;
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

int reconstruction::init(int sizex_sm, int sizey_sm,
                         int sizex_lg, int sizey_lg,
					          int order_sm, int order_lg,
					          const vector<indice>& sten_lg_pre,
					          const vector<double>& mylinwgts_lg,
                         const vector<indice>& sten_sm_pre,
					          const vector<double>& mylinwgts_sm,
                         const MeshInfo& mi, indice start,
					          double mylinwgts_const){

    // clear linear weights
    linwgts_lg.clear();
    linwgts_sm.clear();
    flat_sten_lg.clear();
    flat_sten_sm.clear();

    sten_lg.clear();
    sten_sm.clear();

    stensizelgx = sizex_lg;
    stensizelgy = sizey_lg;

    stensizesmx = sizex_sm;
    stensizesmy = sizey_sm;

    allsize_lgx = mi.MPIglobalCellSize[0] - sizex_lg+1;
    allsize_lgy = mi.MPIglobalCellSize[1] - sizey_lg+1;

    allsize_smx = mi.MPIglobalCellSize[0] - sizex_sm+1;
    allsize_smy = mi.MPIglobalCellSize[1] - sizey_sm+1;

//    for (auto it: sten_lg_pre){
//				indice now = start + it;
//        if (validsten(sizex_lg, sizey_lg, mi, now)){
//				sten_lg.push_back(it);
//				linwgts_lg.push_back(1);
//				flat_sten_lg.push_back(now[1]*allsize_lgx + now[0]);
//        }
//    }

    for (int i=0; i<(int)sten_lg_pre.size(); i++){
        indice now = start + sten_lg_pre.at(i);

        if (validsten(sizex_lg, sizey_lg, mi, now)){

				sten_lg.push_back(sten_lg_pre.at(i));
				linwgts_lg.push_back(mylinwgts_lg.at(i));
				flat_sten_lg.push_back(now[1]*allsize_lgx + now[0]);

		  }
	 }

//    for (auto it: sten_sm_pre){
//				indice now = start + it;
//        if (validsten(sizex_sm, sizey_sm, mi, now)){
//				sten_sm.push_back(it);
//            linwgts_sm.push_back(1);
//				flat_sten_sm.push_back(now[1]*allsize_smx + now[0]);
//        }
//    }

    for (int i=0; i<(int)sten_sm_pre.size(); i++){
        indice now = start + sten_sm_pre.at(i);

        if (validsten(sizex_sm, sizey_sm, mi, now)){

				sten_sm.push_back(sten_sm_pre.at(i));
            linwgts_sm.push_back(mylinwgts_sm.at(i));
				flat_sten_sm.push_back(now[1]*allsize_smx + now[0]);

        }

    }

    if (use_sten_const){
        linwgts_const = mylinwgts_const;
    } else {
        linwgts_const = 0.000;
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

int reconstruction::extractsigma(const vector<double>& sigma_lg, 
                                 const vector<vector<double>>& sigma_sm){

    for (int s=0; s<sten_lg.size(); s++){
        stensigma_lg.at(s) = sigma_lg.at(flat_sten_lg.at(s));
    }

    for (int s=0; s<sten_sm.size(); s++){
        indice position = sten_sm.at(s);
        int po = abs(position[1]) * stensizesmx + abs(position[0]);
        stensigma_sm.at(s) = sigma_sm.at(flat_sten_sm.at(s)).at(po);

//        cout << position[0] << "  " << position[1] << "   "  << sigma_sm.at(s).at(po) <<  endl;
    }

    return 1;
}

int reconstruction::setWgts(double area){

    double sum = 0.0;

    //nonlinwgts_lg.clear(); 
    //nonlinwgts_sm.clear(); 

    for (int l=0; l<nonlinwgts_lg.size(); l++){
        nonlinwgts_lg.at(l) = linwgts_lg.at(l) / pow(abs(stensigma_lg.at(l)) + epsilon*area, s*r_lg + geteta(r_lg));
        sum += nonlinwgts_lg.at(l);
    }

    for (int l=0; l<nonlinwgts_sm.size(); l++){
        nonlinwgts_sm.at(l) = linwgts_sm.at(l) / pow(abs(stensigma_sm.at(l)) + epsilon*area, s*r_sm + geteta(r_sm));
        sum += nonlinwgts_sm.at(l);
    }

    // Initialize the const nonlinear weight
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

    if (!(sum > 0.0) || !std::isfinite(sum)) {
        throw std::runtime_error("Invalid WENO weight normalization.");
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

double reconstruction::eval(double ** localvals, const vertex& p,
                            const vector<tensorstencilpoly>& sten_lg,
                            const vector<tensorstencilpoly>& sten_sm) const{

    double work= 0.0;

    for (int s=0; s<nonlinwgts_lg.size(); s++){
        work += nonlinwgts_lg.at(s) * sten_lg.at(flat_sten_lg.at(s)).eval(localvals,p);
    } 

    for (int s=0; s<nonlinwgts_sm.size(); s++){
        work += nonlinwgts_sm.at(s) * sten_sm.at(flat_sten_sm.at(s)).eval(localvals,p);
    }

    if (use_sten_const){
        work += nonlinwgts_const * localvals[gstart[1]][gstart[0]];
    }

    return work;
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

int reconstruction::printmoreinfo(double ** localval, 
                                  const vector<tensorstencilpoly>& my_sten_lg, 
                                  const vector<tensorstencilpoly>& my_sten_sm) const{

    cout << "Large Stencil coefs : " << endl;
    for (int s=0; s<nonlinwgts_lg.size(); s++){

		  my_sten_lg.at(flat_sten_lg.at(s)).printStencilSol(localval);
        cout << endl;
        my_sten_lg.at(flat_sten_lg.at(s)).printcollapseCoef(localval);
		  cout << endl;
    } 

    cout << endl;

    cout << "Small Stencil coefs : " << endl;
    for (int s=0; s<nonlinwgts_sm.size(); s++){

        my_sten_sm.at(flat_sten_sm.at(s)).printStencilSol(localval);
		  cout << endl;
        my_sten_sm.at(flat_sten_sm.at(s)).printcollapseCoef(localval);
		  cout << endl;
    }


    return 1;
}

int reconstruction::printsigma(){

    cout << "Number of " << sten_lg.size() <<  " large stencil of order : " << r_lg  << " is used." << endl;

    for (int s=0; s<sten_lg.size(); s++ ){
        cout << "At stencil : (" << sten_lg.at(s)[0] + gstart[0] << ", " << 
                                    sten_lg.at(s)[1] + gstart[1] << ") " << 
												stensigma_lg.at(s) << "  " << endl;
    }cout << endl;

    cout << "Number of " << sten_sm.size() <<  " small stencil of order : " << r_sm  << " is used." << endl;

    for (int s=0; s<sten_sm.size(); s++ ){
        cout << "At stencil : (" << sten_sm.at(s)[0] + gstart[0] << ", " << 
                                    sten_sm.at(s)[1] + gstart[1] << ") : " << 
												stensigma_sm.at(s) << "  " << endl;
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
