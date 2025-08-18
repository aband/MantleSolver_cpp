#include "multilevel.h"

// Return whether it is a valid stencil index
// In serial manner
const bool validsten_ml(int sizex, int sizey, 
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

const int geteta_ml(int r){

    if (r==1){
        return 1;
    } else if (r==2){
        return 3;
    } else {
        return 4;
    }
}

// =======================================================================

int singlelevel::init(int order, int sizex, int sizey, const MeshInfo& mi){

    myorder = order;
    mysizex = sizex;
    mysizey = sizey;

    totalsizex = mi.MPIglobalCellSize[0] - sizex +1;
    totalsizey = mi.MPIglobalCellSize[1] - sizey +1;

    totalsten = totalsizex * totalsizey;

    stenpolys.resize(totalsizex * totalsizey);

    for (int j=0; j<totalsizey; j++){
    for (int i=0; i<totalsizex; i++){
        int s = j*totalsizex + i;
        stenpolys.at(s) = tensorstencilpoly(order, sizex, sizey); 
        stenpolys.at(s).setCoef(mi, i, j);
        stenpolys.at(s).setSigma();
        stenpolys.at(s).startx = i;
        stenpolys.at(s).starty = j;
	 }}

    return 1;
}

int mlreconstruction::prepare(const unordered_map<std::string, 
                                                vector<indice>>& methods,
									   const unordered_map<std::string,
										                  indice>& size,
                              indice target, bool use_const, const MeshInfo& mi){

    for (auto& it : methods){
        vector<indice> m;
        my_method.insert(std::make_pair(it.first, m)); 
		  vector<double> lw;
        linwgts.insert(std::make_pair(it.first, lw));
		  vector<double> nlw;
        nonlinwgts.insert(std::make_pair(it.first, nlw));

        for (auto& sten : it.second){
            indice now = target + sten; 

            if (validsten_ml(size.at(it.first)[0], size.at(it.first)[1], mi, now)){
                my_method.at(it.first).push_back(sten);
                linwgts.at(it.first).push_back(1.0);
                nonlinwgts.at(it.first).push_back(1.0);
            }
        }
    }

    if (use_const){
        my_method.insert(std::make_pair<std::string, vector<indice>>
								("const", {{0,0}}));
        linwgts.insert(std::make_pair<std::string, vector<double>>
								("const", {0.001}));
        nonlinwgts.insert(std::make_pair<std::string, vector<double>>
								("const", {0.001}));
    }

    return 1; 
}

int mlreconstruction::printInfo(){

    for (auto& it: my_method){
        cout << "Level : " << it.first << 
					 ". Reconstruciton methods is : " << endl;
		  for (auto& sten: it.second){
            cout << "(" << sten[0] << ", " << sten[1] << "),  ";
        }cout << endl;
    }

    return 1;
}

int mlreconstruction::setWgts(double area){



    return 1;
}

int multilevel::addLevel(int order, const std::string& name,
                         int sizex, int sizey,
                         const MeshInfo& mi){

    // Check no repeating level
    std::unordered_map<std::string, singlelevel>::const_iterator it
				= levels.find(name);
    assert(it == levels.end());

    // Compute number of stencils
    singlelevel level = singlelevel();
    levels.insert(std::make_pair(name, level));
   
    // Compute coefficients and sigam base
    levels.at(name).init(order, sizex, sizey, mi);

    //cout<< name << endl;
    //levels.at(name).stenpolys.at(0).printCoef();
    //cout << endl;

    // Initialization of smoothness indicators
    vector<double> sigma;
    allsigma.insert(std::make_pair(name, sigma));
    allsigma.at(name).resize(levels.at(name).totalsten);
 
    return 1;
}

int multilevel::updateSigma(double ** locvals){

    for (auto& it: levels){

        allsigma.at(it.first).clear();
        allsigma.at(it.first).resize(it.second.totalsten);

        for (int s=0; s<it.second.totalsten; s++){
            allsigma.at(it.first).at(s) = 
						  it.second.stenpolys.at(s).sigma(locvals);
        }

    }

    return 1;
}

int multilevel::printInfo(){

    for (auto& it: levels){
        cout << "At level " << it.first << ": " << endl;
        cout << "This stencil size : " << it.second.mysizex << ", " 
					 << it.second.mysizey << ". Seen as order: "  << it.second.myorder
					 << endl;
		  cout << "Total number of : " << it.second.totalsten << 
					 " has been calculated." << endl;
    }

    return 1;
}
