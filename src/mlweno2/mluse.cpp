#include "mluse.h"

// ================================================================================
int mluse::setmethod(const std::string& pos,
                     const unordered_map<std::string, vector<indice>>& method){

    reconstMethod.erase(pos);
    reconstMethod.insert(std::make_pair(pos, method));

    return 1;
}

int mluse::setbias(const std::string& pos){

    unordered_map<std::string, double> lw;

    for (const auto& it: reconstMethod.at(pos)){
        lw.insert(std::make_pair(it.first, 1));    
    }

    bias.insert(std::make_pair(pos,lw));

    return 1;
}

/*
int mluse::getsol(Tensor<double>& stencilsol, double ** localsol,
                  const indice& stencilindex) const{

    int localcellindexx = 0;
    int localcellindexy = 0;

    for (int j=0; j<stencilsol.getSize(0); j++){
    for (int i=0; i<stencilsol.getSize(1); i++){
        localcellindexx = stencilindex[0] + i;
        localcellindexy = stencilindex[1] + j;
        stencilsol({i,j}) = localsol[localcellindexy][localcellindexx];
    }}

    return 1;
}

int mluse::updatesigma(const multilevel& ml,
                       double ** localsol){

    Tensor<double> stensigma = Tensor<double>(2);
    Tensor<double> stensol   = Tensor<double>(2);
   
    for (const auto& it: ml.reconlevelSet){
        sigma.erase(it);
        int sizex = ml.getSize(it,0); 
        int sizey = ml.getSize(it,1);
        stensigma.setSize({sizex,sizey});
        stensol.setSize({ml.getStencilSize(it,0),ml.getStencilSize(it,1)});
        for (int j=0; j<sizey; j++){
        for (int i=0; i<sizex; i++){
            getsol(stensol, localsol, {i,j}); 
            stensigma({i,j}) = ml.sigma(it, {i,j}, stensol); 
        }}
        sigma.insert(std::make_pair(it, stensigma));
    }

    return 1;
}
*/

int mluse::printsigma(const std::string& name){

    cout << "Print smoothness indicators for reconstruction " << name << endl;

    int sizex = sigma.at(name).getSize(0);
    int sizey = sigma.at(name).getSize(1);

    for (int j=0; j<sizey; j++){
    for (int i=0; i<sizex; i++){
        cout << sigma.at(name)({i,j}) << "  ";
    }cout << endl;}

    return 1;
}

int mluse::geteta(const int& rl) const{
 
    int work = 0;

    if (rl == 1) {
        work = 1;
    } else if (rl == 2) {
        work = 3;
    } else {
        work = 4;
    }

    return work;
}

bool mluse::stencilexist(const multilevel& ml, const indice& index, const std::string& name) const{

    bool work = false;

    if (index[0] >=0 && index[0] < ml.getSize(name, 0) &&
        index[1] >=0 && index[1] < ml.getSize(name, 1)) {

        work = true;
    }

    return work;
}

/*
int mluse::computeWgts(const unordered_map<std::string, vector<indice>>& method,
                       unordered_map<std::string, vector<double>>& wgts,
                       const multilevel& ml, const double& h0,
                       const indice& index){
    indice targetstencilindex;
    double sum = 0.0; 
    int rl = 0;
    int nl = 0;

    for (const auto& it: bias.at("all")){
        // Pick reconstruction levels

        vector<double> nlw;
        nlw.resize(method.at(it.first).size());
        wgts.insert(std::make_pair(it.first,nlw));
        rl = find_max(ml.getStencilSize(it.first,0), ml.getStencilSize(it.first,1));
        nl = geteta(rl);

        for (int m=0; m<method.at(it.first).size(); m++){
            targetstencilindex = index + method.at(it.first).at(m);
            if (stencilexist(ml, targetstencilindex, it.first)) {
                // extract smoothness indicator
                double stensigma = sigma.at(it.first)({targetstencilindex[0],targetstencilindex[1]});
                wgts.at(it.first).at(m) = bias.at("all").at(it.first)/ pow(stensigma + ep * h0*h0, s*rl + nl);
                sum += wgts.at(it.first).at(m);
            } else {
                wgts.at(it.first).at(m) = 0.0;
            }
        }
    }

    for (auto& nw : wgts){
        for (auto& w : nw.second){
            w /= sum;
        }
    }

    return 1;
}
*/

int mluse::computeWgts(const std::string& pos, const multilevel& ml,
                       const indice& index,    const double& h0,
                       unordered_map<std::string, vector<double>>& wgts){

    indice targetstencilindex;
    double sum = 0.0; 
    int rl = 0;
    int nl = 0;

    for (const auto& it: bias.at(pos)){
        // Pick reconstruction levels

        vector<double> nlw;
        nlw.resize(reconstMethod.at(pos).at(it.first).size());
        wgts.insert(std::make_pair(it.first,nlw));
        rl = find_max(ml.getStencilSize(it.first,0), ml.getStencilSize(it.first,1));
        nl = geteta(rl);

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            targetstencilindex = index + reconstMethod.at(pos).at(it.first).at(m);
            if (stencilexist(ml, targetstencilindex, it.first)) {
                // extract smoothness indicator
                double stensigma = ml.getsigma(it.first,{targetstencilindex[0],targetstencilindex[1]});
                wgts.at(it.first).at(m) = bias.at("all").at(it.first)/ pow(stensigma + ep * h0*h0, s*rl + nl);
                sum += wgts.at(it.first).at(m);
            } else {
                wgts.at(it.first).at(m) = 0.0;
            }
        }
    }

    for (auto& nw : wgts){
        for (auto& w : nw.second){
            w /= sum;
        }
    }

    return 1;
}

int mluse::printWgts(const unordered_map<std::string, vector<double>>& wgts){

    for (const auto& it: wgts){
        cout << "Reconstruction level " << it.first << " : ";
        for (const auto& val : it.second){
            cout << val << "  ";
        } cout << endl;
    }

    return 1;
}

/*
double mluse::eval(const vertex& point, const multilevel& ml,
                   const unordered_map<std::string, vector<indice>>& method,
                   const unordered_map<std::string, vector<double>>& wgts,
                   const indice& index, double ** localsol)const {

    double work = 0.0;

    indice targetstencilindex;
 
    for (const auto& it: bias.at("all")){
        // Pick reconstruction levels

        for (int m=0; m<method.at(it.first).size(); m++){
            targetstencilindex = index + method.at(it.first).at(m);
            if (stencilexist(ml, targetstencilindex, it.first)) {
                Tensor<double> sol = Tensor<double>(2);
                sol.setSize({ml.getStencilSize(it.first, 0), ml.getStencilSize(it.first, 1)});
                getsol(sol, localsol, targetstencilindex); 
                // extract smoothness indicator
                work += wgts.at(it.first).at(m) * ml.eval(it.first,{targetstencilindex[0], targetstencilindex[1]},sol,point);
            } else {
                work += 0;
            }
        }
    }

    return work;
}
*/

double mluse::eval(const vertex& point, const multilevel& ml,
                   const std::string& pos, 
                   const unordered_map<std::string, vector<double>>& wgts,
                   const indice& index, double ** localsol) const{

    double work = 0.0;

    indice targetstencilindex;
 
    for (const auto& it: bias.at(pos)){
        // Pick reconstruction levels

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            targetstencilindex = index + reconstMethod.at(pos).at(it.first).at(m);
            if (stencilexist(ml, targetstencilindex, it.first)) {
                Tensor<double> sol = Tensor<double>(2);
                ml.getsol(sol, localsol, targetstencilindex, it.first); 
                // extract smoothness indicator
                work += wgts.at(it.first).at(m) * ml.eval(it.first,{targetstencilindex[0], targetstencilindex[1]},sol,point);
            } else {
                work += 0;
            }
        }
    }

    return work;
}
