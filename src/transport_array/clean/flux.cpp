#include "transport.h"

double LFflux(double fneg, double fpos, double uneg, double upos, double alpha){
    return 0.5*(fneg + fpos - alpha*(upos-uneg)); 
}

int TransportVariable::ExtractThisEdge(const MeshInfo& mi,
                                       indice cellneg, int edgeneg,
                                       indice cellpos, int edgepos,
                                       int gsize,
                                       vector<double>& uneg,
                                       vector<double>& upos){

    uneg.clear();
    upos.clear();

    uneg.resize(gsize);
    upos.resize(gsize);

    // =====================================
    int dofneg = cellneg[1]*mi.MPIglobalCellSize[0] + cellneg[0];
    int dofpos = cellpos[1]*mi.MPIglobalCellSize[0] + cellpos[0];

    for(int g=0; g<gsize; g++){
        uneg.at(g) = my_recon.at(dofneg)->elem_val.at(edgeneg*gsize+g); 
        upos.at(g) = my_recon.at(dofpos)->elem_val.at(edgepos*gsize+g); 
    }

    return 1;
}

static int getNeighbors(int xMaxCell,   int yMaxCell, 
                        int xSize,      int ySize,
                        int i, int j, int& edgepos, int& edgeneg,
                        indice& cellpos, indice& cellneg,
								bool& onbndry){

    // Default setting that the edge locates not on boundary
    onbndry = false;

    if (xSize > xMaxCell){
        // Vertical
        if(i==0){
            // left boundary
            cellneg = {i,j};
            cellpos = {i,j};
            edgepos = 0;
            edgeneg = 0;

            onbndry = true;

        } else if(i==xSize-1){
            // right boundary
            cellneg = {i-1,j};
            cellpos = {i-1,j};
            edgepos = 2;
            edgeneg = 2;

            onbndry = true;

		  } else {
            // interior
            cellneg = {i-1,j};
            cellpos = {i,j};
            edgepos = 0;
            edgeneg = 2;

        }

    } else {
        // Horizontal
        if (j==0){
            // bottom
            cellneg = {i,j};
            cellpos = {i,j};
            edgepos = 1;
            edgeneg = 1;

            onbndry = true;

        } else if (j==ySize-1){
            cellneg = {i,j-1};
            cellpos = {i,j-1};
            edgepos = 3;
            edgeneg = 3;

            onbndry = true;

        } else {
            cellneg = {i,j-1};
            cellpos = {i,j};
            edgepos = 1;
            edgeneg = 3;
        }

    }

    return 1;
}

static int checkEdgeFlux(int xMax, int yMax, const vector<double>& flux){

    for (int j=0; j<yMax; j++){
    for (int i=0; i<xMax; i++){
        cout << "At edge : " << i << ", " << j << " : " << flux.at(j*xMax+i) << "      " ;
    }cout << endl;}

    return 1;
}

static int bndryType(const vertexSet& edge, 
                     const vector<vertex>& vel){

    int type = 0;

    // Assuming velocity semi uniform on each edge
    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);
//cout << unitNormal[0] << "  " << unitNormal[1] << endl;
//cout << vel.at(0)[0] << "  " << vel.at(0)[1] << endl;
    double indicator = unitNormal[0] * vel.at(0)[0] + 
                       unitNormal[1] * vel.at(0)[1];
//cout << indicator << endl;
    if (indicator < -1e-15){
        // inflow boundary
        type = 1;
    } else if (indicator > -1e-15 && indicator < 1e-15){
        // noflow boundary
        type = 2;
    } else if (indicator > 1e-15){
        // outflow boundary
        type = 3;
    }
//cout << type << endl;
    // Type 0 marks an absent boundary type 

    return type; 
}

// advective flux computed at the interior edges
double TransportVariable::advflux_edge(const MeshInfo& mi, 
                                       const vector<vertex>& vel, 
                                       const vector<double>& uneg,
                                       const vector<double>& upos,
                                       const vertexSet& edge,
                                       bool localLF, double gLF){

    double work = 0.0;

    // Compute flux integration along the edge
    const valarray<double>& gwe = GaussWeightsEdge;
	 const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        double fneg = advfunc(uneg.at(g), vel.at(g), unitNormal);
        double fpos = advfunc(upos.at(g), vel.at(g), unitNormal);
//cout <<  vel.at(g)[0] << "  " << vel.at(g)[1] << "  " << unitNormal[0] << "  " << unitNormal[1] << "  " << fneg << "  " << fpos << endl;
        if(localLF){
            gLF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
            gLF = find_max(abs(dfdu(uneg.at(g))),abs(dfdu(upos.at(g))))*gLF;
        }

        work += gwe[g] * LFflux(fpos, fneg, upos.at(g), uneg.at(g), gLF) * len/2.0;
        //work += gwe[g] * LFflux(fpos, fneg, fpos, fneg, gLF) * len/2.0;
    }

    return work;
}

/*
// advective flux computed at the boundary edges
double TransportVariable::advflux_edge(const MeshInfo& mi, 
                                       const vector<vertex>& vel, 
                                       const vector<double>& u,
                                       const vector<double>& f,
                                       const vertexSet& edge,
                                       bool localLF, double gLF){

    double work = 0.0;

    // Compute flux integration along the edge
    const valarray<double>& gwe = GaussWeightsEdge;
	 const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        if(localLF){
            gLF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
            gLF = find_max(abs(dfdu(u.at(g))),abs(dfdu(u.at(g))))*gLF;
        }

        work += gwe[g] * LFflux(f.at(g), f.at(g), u.at(g), u.at(g), gLF) * len/2.0;
    }

    return work;
}
*/

// General function used to compute flux across all the edges
int TransportVariable::advflux_edge_all(const MeshInfo& mi, 
                                        int xSize, int ySize, int offset,
                                        int xMaxCell, int yMaxCell,
                                        const vector<vertex>& edgegaussp,
                                        const vector<vertex>& edgevel,
                                        vector<double>& edgeflux, 
                                        bool localLF, int bndryselect, 
                                        double gLF){

    // Copy edge gauss points and weights
    const valarray<double>& gwe = GaussWeightsEdge;

    int edgepos = 0;
    int edgeneg = 0;

    indice cellpos {0,0};
    indice cellneg {0,0};

    // Velocity across this edge
    vector<vertex> thisedgevel;
    thisedgevel.resize(gwe.size());

    vector<vertex> thisedgegaussp;
    thisedgegaussp.resize(gwe.size());

    // Reconstruction values evaluated from neg cell
    vector<double> uneg;
    // Reconstruction values evaluated from pos cell
    vector<double> upos;
    //vector<double> fneg;
    //vector<double> fpos;

    // Check if this edge is on the boundary or not
    bool onbndry = false;
    int bndrytype = 0;

    // Computed flux value
    double thisflux = 0.0;

    // Loop over all the edges
    for (int j=0; j<ySize; j++){
    for (int i=0; i<xSize; i++){

        getNeighbors(xMaxCell, yMaxCell, xSize, ySize, i, j, 
                     edgepos, edgeneg, cellpos, cellneg, onbndry);

        int dof = j*xSize + i;
        // Edge dof
        dof += offset;

        // Extract velocity from velocity vector
        for (int g=0; g<(int)gwe.size(); g++){
            thisedgevel.at(g) = edgevel.at(dof*gwe.size()+g);
				thisedgegaussp.at(g) = edgegaussp.at(dof*gwe.size()+g); 
        }

		  // Extract four corners vertex of this cell
        vertexSet corners = extractCorners(mi, cellpos);
        int start = (edgepos+3)%4;
        int end   = edgepos;
        // Original edge direction
        vertexSet edge = {corners.at(start), corners.at(end)}; 

        // Uniform edge direction
        vertexSet uniformEdge;

        // Alter boundary vertex order if it is on boundary
        if (onbndry){
            uniformEdge = {edge[1], edge[0]};
        } else {
            uniformEdge = edge;
        }

        // Extract velocity current edge
        ExtractThisEdge(mi, cellneg, edgeneg, 
                            cellpos, edgepos, gwe.size(),
                            uneg, upos);

        // Use uniformEdge to calculate the edge flux
        thisflux = advflux_edge(mi, thisedgevel, uneg, upos, uniformEdge, localLF, gLF);

        //edgeflux.at(dof-offset) = advflux_edge(mi, thisedgevel, 
        //                   uneg, upos, edge, localLF, gLF);

        if (onbndry){

                bndrytype = bndryType(edge, thisedgevel);

                // bndrytype = 1 ------- inflow boundary
                // bndrytype = 2 ------- noflow boundary
					 // bndrytype = 3 ------- outflow boundary

					 if(bndrytype == 1){
                    // If this is not an outflow boundary, the inflow is constrained by dirichlet values
						  // Assign inflow dirichlet values 
                    diriBndry(thisedgegaussp, uneg, bndryselect);

						  thisflux = -1*advflux_edge(mi, thisedgevel, uneg, uneg, uniformEdge, localLF, gLF); 

					 } else if (bndrytype == 2){

                    diriBndry(thisedgegaussp, uneg, bndryselect);

						  thisflux = advflux_edge(mi, thisedgevel, uneg, uneg, uniformEdge, localLF, gLF); 
	
                }
        }

//					 if (xSize > ySize) {cout << "This is a vertical edge" << endl;}
//					 else {cout << "This is a horizontal edge" << endl;}
//                cout << endl;
/*
        if (onbndry){
           cout <<"This edge " << i << ", " <<  j << " is a boundary edge."<< endl;
				cout << "Evaluated flux equals " << thisflux << endl;

            int tmpbt = bndryType(edge, thisedgevel);

				if (tmpbt == 1){
                cout << "This is an inflow boundary" << endl;
				}else if (tmpbt == 2){
                cout << "This is a noflow boundary" << endl;
            } else if (tmpbt == 3){
                cout << "This is a outflow boundary" << endl;
				}

				cout << endl;
        } else {
            //cout <<"This edge " << i << ", " <<  j << " is NOT a boundary edge."<< endl;
        }
*/
        edgeflux.at(dof-offset) = thisflux;

    }}

    return 1;
}

// ======================================================================
/*
double TransportVariable::difflux_edge(){

    double work  = 0.0;

    double len = length

    return work;
}
*/
/*
int TransportVariable::difflux_edge_all(const MeshInfo& mi,
                                        int xSize, int ySize, int offset,
                                        int xMaxCell, int yMaxCell,
                                        vector<double>& edgeflux){
     
    const valarray<double>& gwe = GaussWeightsEdge;

    

    return 1;
}
*/

int TransportVariable::cellflux_all(const MeshInfo& mi, 
                                    double maxv, DM dmu, 
                                    const vector<vertex>& edgegaussp,
                                    const vector<vertex>& edgevel,
												bool globalLF, int bndryselect,
												Vec * influx){

    int N = mi.MPIglobalCellSize[1];
    int M = mi.MPIglobalCellSize[0];

    vector<double> vert;
    vert.resize(N*(M+1));

    vector<double> hori;
    hori.resize(M*(N+1));

    //tolvert = (M+1)*N*3;

    // Vertical flux first
    advflux_edge_all(mi, M+1, N, 0, M, N, edgegaussp, edgevel, vert, globalLF, bndryselect, maxv);

    // horizontal flux next
    advflux_edge_all(mi, M, N+1, (M+1)*N, M, N, edgegaussp, edgevel, hori, globalLF, bndryselect, maxv);

    Vec flux = *influx;

    double ** f;
    PetscCall(DMDAVecGetArray(dmu, flux, &f));

    // Check edge flux
//	 cout << "Horizontal Edge Flux : "  << endl;
//    checkEdgeFlux(M, N+1, hori);

//    cout << endl;

//	 cout << "Vertical Edge Flux : "  << endl;
//    checkEdgeFlux(M+1, N, vert);

//    cout << endl << endl;

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

       double area = mi.cellArea.at(FlatIndic(mi, {i,j}));

       int bottom = j*M+i;
       int top    = (j+1)*M+i;
       int left   = j*(M+1) + i;
       int right  = j*(M+1) + i+1;

       f[j][i] = (hori.at(bottom) - hori.at(top) + vert.at(left) - vert.at(right))/area;

/*
       cout << "This cell " << i << "  " << j << endl
					<<"Bottom : " << hori.at(bottom) << endl 
					<<"Top    : " << hori.at(top) << endl
					<<"Left   : " << vert.at(left) << endl
					<<"Right  : " << vert.at(right)<< endl
					<<"Total  : " << f[j][i] << endl << endl;

    }cout << endl;}
*/
      }}

    DMDAVecRestoreArray(dmu, flux, &f);
 
    return 1;
}
