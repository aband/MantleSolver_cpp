#include "efficient_level.h"

int reconlevel::prepare(const vector<int>& insize,
                        const MeshInfo& mi){
    size = insize;

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

    sp = Tensor<stencilpolynomial>(2);

    sp.setSize({right-left-size[0]+1,
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



    return 1;
}
