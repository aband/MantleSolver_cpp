#include "mluse.h"
#include "util.h"

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

int mluse::setbias(const std::string& pos, 
                   const std::string& level, 
                   const double& newbias){

    bias.at(pos).at(level) = newbias;

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

int mluse::setstencilrange(const MeshInfo& mi){

    int i_start = mi.MPIlocalCellStart[0] - 1;
    int i_end   = mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + 1; 

    int j_start = mi.MPIlocalCellStart[1] - 1;
    int j_end   = mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + 1; 

    int left  = (i_start<0) ? 0 : i_start;
    int right = (i_end > mi.MPIglobalCellSize[0]) ? mi.MPIglobalCellSize[0] : i_end;

    int bottom = (j_start<0) ? 0 : j_start;
    int top    = (j_end > mi.MPIglobalCellSize[1]) ? mi.MPIglobalCellSize[1] : j_end;

    cout << "called ?" << endl << endl;

//    allwgts.setSize({right-left, top-bottom});

    //int converti = 0;
    //int convertj = 0;

    if (left == 0){
        //converti = 0;
        offseti = 0;
    } else {
        //converti = mi.cellGhostLayerSize-1;
        offseti = mi.cellGhostLayerSize-1;
    }

    if (bottom == 0){
        //convertj = 0;
        offsetj = 0;
    } else {
        //convertj = mi.cellGhostLayerSize-1;
        offsetj = mi.cellGhostLayerSize-1;
    }

    return 1;
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
                wgts.at(it.first).at(m) = bias.at(pos).at(it.first)/ pow(stensigma + ep * h0*h0, s*rl + nl);
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

int mluse::computeWgtsConst(const std::string& pos, const multilevel& ml,
                            const indice& index,    const double& h0,
                            weights& wgts){

    indice targetstencilindex;
    double sum = 0;
    for (const auto& it:bias.at(pos)){

        vector<double> nlw;
        nlw.resize(reconstMethod.at(pos).at(it.first).size());
        wgts.insert(std::make_pair(it.first,nlw));

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            targetstencilindex = index + reconstMethod.at(pos).at(it.first).at(m);
            if (stencilexist(ml, targetstencilindex, it.first)) {
                // extract smoothness indicator
                wgts.at(it.first).at(m) = 1;
                sum += 1;
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

// Compute non linear weights with precomputed scaled smoothness indicators
int mluse::computeWgts(const std::string& pos, const indice& gcell, const multilevel& ml, weights& wgts){

    indice targetstencilindex;
    double sum = 0;

    for (const auto& it: bias.at(pos)){
        vector<double> nlw;
        nlw.resize(reconstMethod.at(pos).at(it.first).size());
        wgts.insert(std::make_pair(it.first,nlw));
 
        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            targetstencilindex = gcell + reconstMethod.at(pos).at(it.first).at(m);
            if (stencilexist(ml, targetstencilindex, it.first)) {
                // extract smoothness indicator
                wgts.at(it.first).at(m) = bias.at(pos).at(it.first)* ml.getscaledsigma(it.first,{targetstencilindex[0],
                                                                                                 targetstencilindex[1]});

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

int mluse::computeWgtsConst(const multilevel& ml, 
                            const MeshInfo& mi,
                            const double& h0,
                            Tensor<weights>& allwgts){

    // give a constant averaged nonlinear weights
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

            computeWgtsConst(pos(mi,stencilindex), ml, stencilindex, h0, allwgts({i,j}));
        }
    }


    return 1;
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

int mluse::computeWgts(const multilevel& ml, const MeshInfo& mi, Tensor<weights>& allwgts){

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

            computeWgts(pos(mi,stencilindex), stencilindex, ml, allwgts({i,j}));
        }
    }


    return 1;
}

// Using different combination for different positions
int mluse::computeWgts(const multilevel& ml, 
                       const MeshInfo& mi, 
							  const double& h0,
                       Tensor<weights>& allwgts,
                       posFunc pfunc){

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

            computeWgts(pfunc(mi,{i,j}), ml, stencilindex, h0, allwgts({i,j}));
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
                // extract nonlinear weights
                work += wgts.at(it.first).at(m) * ml.eval(it.first,{targetstencilindex[0], targetstencilindex[1]},sol,point);
            } else {
                work += 0;
            }
        }
    }

    return work;
}

int mluse::eff_order(const indice& index,    const multilevel& ml, 
                     const std::string& pos, const weights& wgts) const{

    int work = 0;

    indice targetstencilindex;

    for (const auto& it: bias.at(pos)){

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            targetstencilindex = index + reconstMethod.at(pos).at(it.first).at(m);
            if (stencilexist(ml, targetstencilindex, it.first)){

                work += wgts.at(it.first).at(m) * ml.getorder(it.first);
					 cout << ml.getorder(it.first) << endl;
            }
        }
    }

    return work;
}

int mluse::sumscaled(const multilevel& ml,
                     double& sum, 
                     derivative& sumder,
                     const std::string& pos, 
                     const indice& gcell) const{

    indice targetstencilindex;
    sum = 0.0;
    sumder.clear();

    for (const auto& it: bias.at(pos)){

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            targetstencilindex = gcell + reconstMethod.at(pos).at(it.first).at(m);

            if (stencilexist(ml, targetstencilindex, it.first)) {

                sum += bias.at(pos).at(it.first)* 
                       ml.getscaledsigma(it.first, {targetstencilindex[0], targetstencilindex[1]});
                unordered_map_arithmetic(sumder, ml.getderscaledsigma(it.first,{targetstencilindex[0],targetstencilindex[1]}), 
                                         std::plus<double>(), bias.at(pos).at(it.first), std::multiplies<double>());
            }
        }
    }

    return 1;
}

/*!
 * Stencil index equal to local index of its left bottom cell.
 */

int mluse::der(const vertex& point, const multilevel& ml,
               const::string& pos,  
               const indice& index, double ** localsol,
               const MeshInfo& mi,
               unordered_map<int, double>& der) const{

    // Derivative of a reconstruction consists of two parts
    // dR/du = \sum dw/du P = \sum w dP/du

    // Clear target derivative object
    der.clear();

    indice targetstencilindex;

    double sum;
    derivative sumder;

    // Compute sum of weights and sum of derivative of weights
    sumscaled(ml,sum, sumder, pos, index);

    for (const auto& it : bias.at(pos)){
        // Loop through all levels first

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            // Loop through method to find target stencil

            targetstencilindex = index + reconstMethod.at(pos).at(it.first).at(m);

            if (stencilexist(ml, targetstencilindex, it.first)) {
                // Check if this stencil actually exists

                Tensor<double> sol = Tensor<double>(2);
                ml.getsol(sol, localsol, targetstencilindex, it.first); 

                // Extract non linear weight for this tencil
                double scaled = bias.at(pos).at(it.first)*
                ml.getscaledsigma(it.first, {targetstencilindex[0], 
                                             targetstencilindex[1]});

                double nlw = scaled/sum;

                // Compute dw/du for each stencil
                derivative dscaled = ml.getderscaledsigma(it.first,
                {targetstencilindex[0], targetstencilindex[1]});
/*
                cout << "dsigmas : " << endl;
                unordered_map_print(ml.getdersigma(it.first, {targetstencilindex[0], targetstencilindex[1]}));

                cout << "scaled sigmas : " << endl;
                cout << ml.getscaledsigma(it.first, {targetstencilindex[0], targetstencilindex[1]}) << endl;

                cout << "sigma : " << endl;
                cout << ml.getsigma(it.first, {targetstencilindex[0], targetstencilindex[1]}) << endl;

                cout << "dscalde sigmas : " << endl;
                unordered_map_print(ml.getderscaledsigma(it.first, {targetstencilindex[0], targetstencilindex[1]}));
*/
                unordered_map_arithmetic(dscaled, bias.at(pos).at(it.first),
                                         std::multiplies<double>());

                unordered_map_arithmetic(dscaled, 1.0/sum, std::multiplies<double>()); 
                unordered_map_arithmetic(dscaled, sumder, std::plus<double>(), 
                -1*scaled/sum/sum, std::multiplies<double>());

                double val = ml.eval(it.first, 
                {targetstencilindex[0],targetstencilindex[1]}, sol, point);

                // p1 = dw/du * p
                unordered_map_arithmetic(dscaled, val, std::multiplies<double>()); 

                for (int j=0; j<sol.getSize(1); j++){
                for (int i=0; i<sol.getSize(0); i++){
                    // Loop through the selected solution stencil 
                    // Be careful with parallel there will be offset for paralle case
                    indice globalcell = {i,j}; // Reterive global cell index first
                    globalcell += targetstencilindex;
                    int flatgcell = FlatIndic(mi, globalcell);

                    // The complete derivative consists of two parts
                    // p representing dp/du
                    double p = 0.0;

                    // p = w * dp/du
                    p = nlw * ml.eval(it.first, {targetstencilindex[0], targetstencilindex[1]},
                                                 {i,j}, point);

                    std::unordered_map<int, double>::const_iterator got = dscaled.find(flatgcell);

                    if (got == dscaled.end()){
                        // This derivative has not been calculated
                        dscaled.insert(std::make_pair(flatgcell, p));
                    } else {
                        dscaled.at(flatgcell) += p;
                    }

                    //der[flatgcell] += p;
                }}

                unordered_map_arithmetic(der,dscaled,std::plus<double>());
            } 
        }
    }
    return  1;
}

// For testing purpose
int mluse::derpseudo(const vertex& point, const multilevel& ml,
                     const std::string& pos, const weights& wgts,
                     const indice& index, double ** localsol,
                     const MeshInfo& mi,
                     derivative& der) const{

    // Clear target derivative object
    der.clear();

    indice targetstencilindex;

    for (const auto& it: bias.at(pos)){

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){

            targetstencilindex = index + reconstMethod.at(pos).at(it.first).at(m);

            if (stencilexist(ml, targetstencilindex, it.first)) {
 
                Tensor<double> sol = Tensor<double>(2);
                ml.getsol(sol, localsol, targetstencilindex, it.first); 

                double nlw = wgts.at(it.first).at(m);

                for (int j=0; j<sol.getSize(1); j++){
                for (int i=0; i<sol.getSize(0); i++){
                    // Loop through the selected solution stencil 
                    // Be careful with parallel there will be offset for paralle case
                    indice globalcell = {i,j}; // Reterive global cell index first
                    globalcell += targetstencilindex;
                    int flatgcell = FlatIndic(mi, globalcell);

                    // The complete derivative consists of two parts
                    // p representing dp/du
                    double p = 0.0;

                    // p = w * dp/du
                    p = nlw * ml.eval(it.first, {targetstencilindex[0], targetstencilindex[1]},
                                                 {i,j}, point);

                    std::unordered_map<int, double>::const_iterator got = der.find(flatgcell);

                    if (got == der.end()){
                        // This derivative has not been calculated
                        der.insert(std::make_pair(flatgcell, p));
                    } else {
                        der.at(flatgcell) += p;
                    }
                }}
            }
        }
    }

    return 1;
}

int mluse::dnlwtest(const multilevel& ml,
                    const std::string& pos,
                    const indice& index, 
                    double ** localsol,
                    const MeshInfo& mi) const{

    // Derivative of a reconstruction consists of two parts
    // dR/du = \sum dw/du P = \sum w dP/du

    indice targetstencilindex;

    double sum;
    derivative sumder;

    // Compute sum of weights and sum of derivative of weights
    sumscaled(ml,sum, sumder, pos, index);

    for (const auto& it : bias.at(pos)){
        // Loop through all levels first

        for (int m=0; m<reconstMethod.at(pos).at(it.first).size(); m++){
            // Loop through method to find target stencil

            targetstencilindex = index + reconstMethod.at(pos).at(it.first).at(m);

            if (stencilexist(ml, targetstencilindex, it.first)) {
                // Check if this stencil actually exists

                Tensor<double> sol = Tensor<double>(2);
                ml.getsol(sol, localsol, targetstencilindex, it.first); 

                // Extract non linear weight for this tencil
                double scaled = bias.at(pos).at(it.first)*
                ml.getscaledsigma(it.first, {targetstencilindex[0], 
                                             targetstencilindex[1]});

                double nlw = scaled/sum;

                // Compute dw/du for each stencil
                derivative dscaled = ml.getderscaledsigma(it.first,
                {targetstencilindex[0], targetstencilindex[1]});

                unordered_map_arithmetic(dscaled, bias.at(pos).at(it.first),
                                         std::multiplies<double>());

                unordered_map_arithmetic(dscaled, 1.0/sum, std::multiplies<double>()); 
					 cout << endl;
unordered_map_print(dscaled);

                unordered_map_arithmetic(dscaled, sumder, std::plus<double>(), 
                -1*scaled/sum/sum, std::multiplies<double>());

derivative testsumder = sumder;
                unordered_map_arithmetic(testsumder, -1*scaled/sum/sum, std::multiplies<double>());
unordered_map_print(testsumder);

                cout << it.first << " stencil index : " << targetstencilindex[0] << "  " << targetstencilindex[1] << endl;
                unordered_map_print(dscaled);

            } 
        }
    }

    return 1;
}
