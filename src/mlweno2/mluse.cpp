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

int mluse::computeWgts(const std::string& pos, const multilevel& ml,
                       const indice& index,    const double& h0,
                       weights& wgts){

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

int mluse::printWgts(const weights& wgts){

    for (const auto& it: wgts){
        cout << "Reconstruction level " << it.first << " : ";
        for (const auto& val : it.second){
            cout << val << "  ";
        } cout << endl;
    }

    return 1;
}

std::string pos(const MeshInfo& mi, const indice& stencilindex){

    std::string position = "all"; // For testing

    return position;
}

int mluse::computeWgts(const multilevel& ml, const MeshInfo& mi,const double& h0, 
                       Tensor<weights>& allwgts){

    allwgts = Tensor<weights>(2);

    int i_start = mi.MPIlocalCellStart[0] - 1;
    int i_end   = mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + 1; 

    int j_start = mi.MPIlocalCellStart[1] - 1;
    int j_end   = mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + 1; 

    int left  = (i_start<0) ? 0 : i_start;
    int right = (i_end > mi.MPIglobalCellSize[0]) ? mi.MPIglobalCellSize[0] : i_end;

    int bottom = (j_start<0) ? 0 : j_start;
    int top    = (j_end > mi.MPIglobalCellSize[1]) ? mi.MPIglobalCellSize[1] : j_end;

    allwgts.setSize({right-left, top-bottom});

    int converti = 0;
    int convertj = 0;

    if (left == 0){
        converti = 0;
    } else {
        converti = mi.cellGhostLayerSize-1;;
    }

    if (bottom == 0){
        convertj = 0;
    } else {
        convertj = mi.cellGhostLayerSize-1;;
    }

    for (int j=0; j<top-bottom; j++){
        for (int i=0; i<right-left; i++){
            weights wgts;
            // convert index
            indice stencilindex {i + converti, 
                                 j + convertj};

            computeWgts(pos(mi,stencilindex), ml, stencilindex, h0, allwgts({i,j}));
        }
    }

    return 1;
}

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
