#include "reconstruction.h"

int reconstruction::prepare(const vector<int>& insize, const MeshInfo& mi){

    size = insize;

    stencilPoly = Tensor<stencilpolynomial>(2);

    // Determine ranges
    int left, right, top, bottom;

    int j_start = mi.MPIlocalCellStart[1] - mi.cellGhostLayerSize;
    int i_start = mi.MPIlocalCellStart[0] - mi.cellGhostLayerSize;
    int j_end   = mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + 
                  mi.cellGhostLayerSize; 
    int i_end   = mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + 
                  mi.cellGhostLayerSize; 

    left  = (i_start<0) ? 0 : i_start;
    right = (i_end > mi.MPIglobalCellSize[0]) ? mi.MPIglobalCellSize[0] : i_end;

    bottom = (j_start<0) ? 0 : j_start;
    top    = (j_end > mi.MPIglobalCellSize[1]) ? mi.MPIglobalCellSize[1] : j_end;

    stencilPoly = Tensor<stencilpolynomial>(2);

    stencilPoly.setSize({right-left-size[0]+1,
                         top-bottom-size[1]+1});

    // Using average values for this scale
	 // and area value
    double area  = mi.L*mi.H/(double)(mi.MPIglobalCellSize[0] *
                                      mi.MPIglobalCellSize[1]);
    double scale = sqrt(area);

    // Initialize stencil polynomials 
    // Compute coefficients and its corresponding sigmas
    for (int s=0; s<stencilPoly.getSize(1); s++){
    for (int k=0; k<stencilPoly.getSize(0); k++){
        // Extract vector of four corners of the cells in the stencil
        vector<vector<vertex>> cornerSet; 
        vector<vertex> refcell;
        vertex center;
        for (int j=0; j<size[1]; j++){
            for (int i=0; i<size[0]; i++){
                indice global {i+mi.MPIlocalCellStart[0] + k, 
                               j+mi.MPIlocalCellStart[1] + s};
                vector<vertex> cellCornerSet = extractCorners(mi, global);
                cornerSet.push_back(cellCornerSet);
            } 
        }

        getcenter(cornerSet, refcell, center, scale, area);
        stencilPoly({k,s}) = stencilpolynomial(size[0],size[1]);
        stencilPoly({k,s}).center = center;
        stencilPoly({k,s}).h = scale;
        stencilPoly({k,s}).setCoef(cornerSet, center, scale);
        //stencilPoly({k,s}).printCoef();
        stencilPoly({k,s}).sigma(refcell, area, center, scale);
    }}

    // Test print output ==========================================================
/*
    cout << "Test stencil size is " << size[0] << " " << size[1] << endl;
    cout << "Untrimmed domain : " << endl 
         << " x in " << i_start << " : "  << i_end << endl
         << " y in " << j_start << " : "  << j_end << endl;
    cout << "Stencils are created in the trimmed domain : " << endl 
         << " x in " << left << " : "  << right << endl
         << " y in " << bottom << " : "  << top << endl;

    cout << "There are " << stencilPoly.getSize(0) << " x " 
                         << stencilPoly.getSize(1) << " stencils created." << endl;
*/
    // ============================================================================
    return 1;
}

// Get center point , reference cell corners and scale factors 
int reconstruction::getcenter(const vector<vector<vertex>>& cornerSet,
                              vector<vertex>& refcell,
                              vertex& center, double& h, double& area){

    Tensor<int> tmp = Tensor<int>(2);
    tmp.setSize(size);

    center = (cornerSet.at(0)[0] + cornerSet.at(size[0]-1)[1] +
              cornerSet.at(size[0]*size[1]-1)[2] + 
              cornerSet.at(tmp.getIndex({0,size[1]-1}))[3] )/4.0;

    // Create reference cell corners
    refcell.resize(4);
    vertex add = {-h/2, -h/2};
    refcell[0] = center + add;
    add = {h/2, -h/2};
    refcell[1] = center + add;
    add = {h/2, h/2};
    refcell[2] = center + add;
    add = {-h/2, h/2};
    refcell[3] = center + add; 

    return 1;
}

double reconstruction::eval(const vector<int>& index, 
                            const Tensor<double>& sol,
                            const vertex& point) const{
    return stencilPoly(index).eval(sol, point);
}

double reconstruction::sigma(const vector<int>& index,
                             const Tensor<double>& sol) const{

    return stencilPoly(index).sigma(sol);
}

// ================================================================================
int multilevel::addLevel(const std::string& name, 
                         const vector<int>& stencilSize,
                         const MeshInfo& mi){

    std::unordered_map<std::string, reconstruction>::const_iterator it = 
    mlrecons.find(name);
    // Check if the level existing
    assert(it == mlrecons.end());

    reconstruction level = reconstruction();    
    level.prepare(stencilSize, mi); 

    mlrecons.insert(std::make_pair(name, level));

    reconlevelSet.insert(name);

    return 1;
}

double multilevel::eval(const std::string& name, 
                        const vector<int>& index, 
                        const Tensor<double>& sol,
                        const vertex& point) const{

    return mlrecons.at(name).eval(index,sol,point);
} 

double multilevel::sigma(const std::string& name, 
                         const vector<int>& index,
                         const Tensor<double>& sol) const{

    return mlrecons.at(name).sigma(index, sol);
}

int multilevel::getsol(Tensor<double>& stencilsol, double ** localsol,
                       const indice& stencilindex, const std::string& name) const{

    int localcellindexx = 0;
    int localcellindexy = 0;

    stencilsol.setSize({getStencilSize(name, 0), getStencilSize(name, 1)});

    for (int j=0; j<stencilsol.getSize(0); j++){
    for (int i=0; i<stencilsol.getSize(1); i++){
        localcellindexx = stencilindex[0] + i;
        localcellindexy = stencilindex[1] + j;
        stencilsol({i,j}) = localsol[localcellindexy][localcellindexx];
    }}

    return 1;
}

int multilevel::updatesigma(double ** localsol){

    Tensor<double> stensigma = Tensor<double>(2);
    Tensor<double> stensol   = Tensor<double>(2);
   
    for (const auto& it: reconlevelSet){
        alllevelsigma.erase(it);
        int sizex = getSize(it,0); 
        int sizey = getSize(it,1);
        stensigma.setSize({sizex,sizey});
        //stensol.setSize({getStencilSize(it,0),getStencilSize(it,1)});
        for (int j=0; j<sizey; j++){
        for (int i=0; i<sizex; i++){
            getsol(stensol, localsol, {i,j}, it); 
            stensigma({i,j}) = sigma(it, {i,j}, stensol); 
        }}
        alllevelsigma.insert(std::make_pair(it, stensigma));
    }

    return 1;
}

int multilevel::printsigma(const std::string& name){

    cout << "Print smoothness indicators for reconstruction " << name << endl;

    int sizex = alllevelsigma.at(name).getSize(0);
    int sizey = alllevelsigma.at(name).getSize(1);

    for (int j=0; j<sizey; j++){
    for (int i=0; i<sizex; i++){
        cout << alllevelsigma.at(name)({i,j}) << "  ";
    }cout << endl;}

    return 1;
}

int multilevel::printcoef(const std::string& name){

    cout << "Print stencil polynomial coefficients for reconstruction " << name << endl;

    mlrecons.at(name).printcoef();

    return 1;
}
