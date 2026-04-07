#include "transport.h"

// Create Default (3,2) reconstruction
int TransportVariable::CreateDefaultReconstruction(const MeshInfo& mi){


	 int sizelgx = 3;
    int sizelgy = 3;
    int orderlg = 2;

	 int sizesmx = 2;
    int sizesmy = 2;
    int ordersm = 1;

    // (3,2) reconstruction but 1D
    vector<indice> sten_lg_pre = {{-1,-1}};
    vector<indice> sten_sm_pre = {{0,-1}, {0,0}, {-1,-1}, {-1,0}};


/*
	 int sizelgx = 5;
    int sizelgy = 5;
    int orderlg = 4;

	 int sizesmx = 3;
    int sizesmy = 3;
    int ordersm = 2;

    // (3,2) reconstruction but 1D
    vector<indice> sten_lg_pre = {{-2,-2}};
    vector<indice> sten_sm_pre = {{0,-2}, {0,0}, {-2,-2}, {-2,0}};
*/

    CreateReconstruction(mi, sizelgx, sizelgy, orderlg, 
                             sizesmx, sizesmy, ordersm,
                             sten_lg_pre, sten_sm_pre,
									  false); 

    return 1;
}

int TransportVariable::CreateReconstruction(const MeshInfo& mi, 
                                            int sizelgx, int sizelgy, int orderlg,
                                            int sizesmx, int sizesmy, int ordersm,
                                            vector<indice>& sten_lg_pre,
                                            vector<indice>& sten_sm_pre,
														  bool use_sten_const){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    int Mlg = M-sizelgx+1;
    int Nlg = N-sizelgy+1;

    int Msm = M-sizesmx+1;
    int Nsm = N-sizesmy+1;

    stenlg.resize(Mlg*Nlg);

    for (int j=0; j<Nlg; j++){
    for (int i=0; i<Mlg; i++){
        int s = j*Mlg+i;
        stenlg.at(s) = tensorstencilpoly(orderlg, sizelgx, sizelgy);
        stenlg.at(s).setCoef(mi,i,j);
		  stenlg.at(s).setSigma();
		  stenlg.at(s).startx = i;
		  stenlg.at(s).starty = j;
    }}

    stensm.resize(Msm*Nsm);

    for (int j=0; j<Nsm; j++){
    for (int i=0; i<Msm; i++){
        int s = j*Msm + i;
        stensm.at(s) = tensorstencilpoly(ordersm, sizesmx, sizesmy);
        stensm.at(s).setCoef(mi,i,j);
		  stensm.at(s).setSigma();
		  stensm.at(s).startx = i;
		  stensm.at(s).starty = j;
    }}

    my_recon.resize(M*N);

    // Initializing reconstrucitons for each cell
	 // Change here to have sided reconstruction on the boundary
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
 
        my_recon.at(s) = new reconstruction();

		  my_recon.at(s)->use_sten_const = use_sten_const;

        my_recon.at(s)->init(sizesmx,sizesmy,
                             sizelgx,sizelgy,
                             ordersm,orderlg,
                             sten_lg_pre, sten_sm_pre, mi,{i,j});
    }}

    // Allocate memory for extra evaluations
    cellgauss.resize(M*N*9);
    cellcenter.resize(M*N);

    // Default diffusion is turned off
    if (diffusion == true){
//        cout << "Diffusion is defined" << endl;
        //degree = 3 + 1;
        degree = 2 + 1;
        halfpts = std::ceil((degree+1)/2.0);
        numpts  = halfpts * 2;

        lagDer.init(numpts - 1);

        samplingp.resize(((mi.MPIglobalCellSize[0]+1)*mi.MPIglobalCellSize[1]+
                           mi.MPIglobalCellSize[0]*(mi.MPIglobalCellSize[1]+1))*numpts);

        samplingv.resize(((mi.MPIglobalCellSize[0]+1)*mi.MPIglobalCellSize[1]+
                           mi.MPIglobalCellSize[0]*(mi.MPIglobalCellSize[1]+1))*numpts);

	 }

    return 1;
}

int TransportVariable::CreateReconstruction(const MeshInfo& mi, 
                                            int sizelgx, int sizelgy, int orderlg,
                                            int sizesmx, int sizesmy, int ordersm,
                                            vector<indice>& sten_lg_pre,
														  vector<double>& mylinwgts_lg,
                                            vector<indice>& sten_sm_pre,
														  vector<double>& mylinwgts_sm,
														  bool use_sten_const,
														  double mylinwgts_const){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    int Mlg = M-sizelgx+1;
    int Nlg = N-sizelgy+1;

    int Msm = M-sizesmx+1;
    int Nsm = N-sizesmy+1;

    stenlg.resize(Mlg*Nlg);

    for (int j=0; j<Nlg; j++){
    for (int i=0; i<Mlg; i++){
        int s = j*Mlg+i;
        stenlg.at(s) = tensorstencilpoly(orderlg, sizelgx, sizelgy);
        stenlg.at(s).setCoef(mi,i,j);
		  stenlg.at(s).setSigma();
		  stenlg.at(s).startx = i;
		  stenlg.at(s).starty = j;
    }}

    stensm.resize(Msm*Nsm);

    for (int j=0; j<Nsm; j++){
    for (int i=0; i<Msm; i++){
        int s = j*Msm + i;
        stensm.at(s) = tensorstencilpoly(ordersm, sizesmx, sizesmy);
        stensm.at(s).setCoef(mi,i,j);
		  stensm.at(s).setSigma();
		  stensm.at(s).startx = i;
		  stensm.at(s).starty = j;
    }}

    my_recon.resize(M*N);

    // Initializing reconstrucitons for each cell
	 // Change here to have sided reconstruction on the boundary
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
 
        my_recon.at(s) = new reconstruction();

		  my_recon.at(s)->use_sten_const = use_sten_const;

        my_recon.at(s)->init(sizesmx,sizesmy,
                             sizelgx,sizelgy,
                             ordersm,orderlg,
                             sten_lg_pre, mylinwgts_lg,
									  sten_sm_pre, mylinwgts_sm, 
									  mi,{i,j}, mylinwgts_const);
    }}

    // Allocate memory for extra evaluations
    cellgauss.resize(M*N*9);
    cellcenter.resize(M*N);

    // Default diffusion is turned off
    if (diffusion == true){
//        cout << "Diffusion is defined" << endl;
        //degree = 3 + 1;
        degree = 3 + 1;
        halfpts = std::ceil((degree+1)/2.0);
        numpts  = halfpts * 2;

        lagDer.init(numpts - 1);

        samplingp.resize(((mi.MPIglobalCellSize[0]+1)*mi.MPIglobalCellSize[1]+
                           mi.MPIglobalCellSize[0]*(mi.MPIglobalCellSize[1]+1))*numpts);

        samplingv.resize(((mi.MPIglobalCellSize[0]+1)*mi.MPIglobalCellSize[1]+
                           mi.MPIglobalCellSize[0]*(mi.MPIglobalCellSize[1]+1))*numpts);

	 }

    return 1;
}

// Update reconstruction nonlinear weights and smoothness indicators
int TransportVariable::UpdateRecon(const MeshInfo& mi, double ** locvals){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    sigma_lg.clear();
    sigma_lg.resize(stenlg.size());

    sigma_sm.clear();
    sigma_sm.resize(stensm.size());

    for (int s=0; s<stenlg.size(); s++){
        sigma_lg.at(s) = stenlg.at(s).sigma(locvals);
    }

    for (int s=0; s<stensm.size(); s++){
       sigma_sm.at(s) = stensm.at(s).sigma(locvals);
    }

    // Setup nonlinear weights
    for (int s=0; s<my_recon.size(); s++){
        my_recon.at(s)->extractsigma(sigma_lg, sigma_sm);
        my_recon.at(s)->setWgts(1.0/(double)M/(double)N);
    }

    return 1;
}

// Evaluate at each gauss points on the edges
int TransportVariable::Evaluate(const MeshInfo& mi, DM dmu){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    Vec localu;
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, sol, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, sol, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    UpdateRecon(mi, lu);

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        //EvaluateEdge(i, j, mi, lu);
        EvaluateExtra(i, j, mi, lu);

    }}

    if (diffusion == true){
       EvaluateSamplesEdge(mi, M+1, N  , 0      , M, N, lu);
       EvaluateSamplesEdge(mi, M  , N+1, (M+1)*N, M, N, lu);
	 }

    PrintPatch(mi, 10, 10, lu, "patch.dat", "patchx.dat", "patchy.dat");

    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}

int TransportVariable::getNeighbors(int xMaxCell,   int yMaxCell, 
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

int TransportVariable::getUniformEdge(const vertexSet& edge, vertexSet& uniformEdge, const indice& cellid,
					                       int xSize, int ySize, int xMaxCell, int yMaxCell){

    if (xSize > xMaxCell){
        // Vertical
        if (cellid[0] == 0){
            uniformEdge = edge;
        } else {
            uniformEdge = {edge[1], edge[0]};
		  }
    } else {
        // Horizontal
        if (cellid[1] == 0){
            uniformEdge = edge;
		  } else {
            uniformEdge = {edge[1], edge[0]};
		  }
	 }

//    if (cellid[0] == 0 || cellid[1] == 0){
//	 } else {
//   }

    return 1;
}

// Evaluate reconstruction value at each gauss points on the edges
// With the same order left - bottom - right - top
int TransportVariable::EvaluateEdge(int i, int j, const MeshInfo& mi, double ** locvals){

    // Get gauss points on the edges 
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    vertexSet corners = extractCorners(mi, {i,j});

    vector<vertex> gaussp;
    gaussp.resize(4*gpe.size());

    vector<indice> edgeBound {{3,0},{0,1},{2,1},{3,2}};

    for (int e=0; e<4; e++){
        vertexSet edge = {corners.at(edgeBound.at(e)[0]), corners.at(edgeBound.at(e)[1])};

        for (int g=0; g<gpe.size(); g++){
            gaussp.at(e*gpe.size()+g) = GaussMapPointsEdge({gpe[g]}, edge);
        }
    }

    my_recon.at(j*mi.MPIglobalCellSize[0]+i)->eval(locvals, gaussp, stenlg, stensm); 

    return 1;
}

// Evaluate gauss points on the edges and inside the cells
int TransportVariable::EvaluateExtra(int i, int j, const MeshInfo& mi, double ** locvals){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    vertexSet corners = extractCorners(mi, {i,j});

    vector<vertex> gaussp;
    gaussp.resize(4*gpe.size());

    // Evaluate at edge points
    vector<indice> edgeBound {{3,0},{0,1},{2,1},{3,2}};

    int cellID = j*mi.MPIglobalCellSize[0]+i;

    for (int e=0; e<4; e++){
        vertexSet edge = {corners.at(edgeBound.at(e)[0]), corners.at(edgeBound.at(e)[1])};

        for (int g=0; g<gpe.size(); g++){
            gaussp.at(e*gpe.size()+g) = GaussMapPointsEdge({gpe[g]}, edge);
        }
    }

    my_recon.at(cellID)->eval(locvals, gaussp, stenlg, stensm); 

    // Evaluate at cell points
    for (int g=0; g<(int)gwf.size(); g++){
        vertex cellmapped = GaussMapPointsFace(gpf[g], corners);
        cellgauss.at(cellID*(int)gwf.size() + g) = my_recon.at(cellID)->eval(locvals,cellmapped,stenlg,stensm);
    }	

    // Evaluate at cell center
    vertex cellcentergrid = GaussMapPointsFace({0.0,0.0}, corners);
    cellcenter.at(cellID) = my_recon.at(cellID)->eval(locvals,cellcentergrid,stenlg,stensm);

    return 1;
}

// Extract sampling values
int TransportVariable::EvaluateSamples(const vertexSet& edge,
                                       const vertex& unitNormal,
                                       double dx, double len,
										 		   int halfpts,
										 		   int cellidneg, int cellidpos,
										 		   double ** locvals,
                                       vector<double>& samples,
													vector<vertex>& samplesv){

//    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    samples.clear();
	 samples.resize((int)gpe.size() * numpts);

    samplesv.clear();
	 samplesv.resize((int)gpe.size() * numpts);

    for (int g=0; g<(int) gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        for (int s=0; s<halfpts; s++){

            //vertex point0 = mapped - (halfpts - 0.5 - s)*dx*unitNormal; 
            //vertex point1 = mapped + (0.5+s)*dx*unitNormal;
            vertex point0 = mapped - (halfpts - 0.5 - s)*dx*unitNormal; 
            vertex point1 = mapped + (0.5+s)*dx*unitNormal;

            samples.at(g*2*halfpts + s) =  
						  my_recon.at(cellidneg)->eval(locvals, point0, stenlg, stensm)/dx*len/2.0;
            samples.at(g*2*halfpts + halfpts + s) = 
						  my_recon.at(cellidpos)->eval(locvals, point1, stenlg, stensm)/dx*len/2.0;

            samplesv.at(g*2*halfpts + s) = point0; 
            samplesv.at(g*2*halfpts + halfpts + s) = point1; 

		  }

	 }

    return 1;
}

// Evaluate samples but for edges on boundary
int TransportVariable::EvaluateSamples(const MeshInfo& mi, 
                                       const indice& cellid,
													const vertexSet& edge,
													const vertex& unitNormal,
													double dx, double len,
													int xSize, int ySize, 
													int xMaxCell, int yMaxCell,
													double ** locvals,
													vector<double>& samples,
													vector<vertex>& samplesv){

    const valarray<double>& gpe = GaussPointsEdge;
    samples.clear();
	 samples.resize((int)gpe.size() * numpts);

    samplesv.clear();
	 samplesv.resize((int)gpe.size() * numpts);

    for (int g=0; g<(int)gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);


        if (diffBndryType(mapped) == 1){

            for (int s=0; s<halfpts; s++){
                if (xSize > xMaxCell){
                    // Vertical
                    if (cellid[0] == 0){
                        // Left sampling values should be replaced with dirichlet values
                        vertex point0 = mapped - (halfpts - 0.5 - s)*dx*unitNormal; 
                        vertex point1 = mapped + (0.5+(s-1))*dx*unitNormal;
    
                        samples.at(g*2*halfpts + s) = 
		    			   	 my_recon.at(FlatIndic(mi,cellid))->eval(locvals, point0, stenlg, stensm)/dx*len/2.0;
//                        samples.at(g*2*halfpts + halfpts + s) = my_recon.at(cellid)->eval(locvals, point1, stenlg, stensm);
                        samples.at(g*2*halfpts + halfpts + s) = diffBndry(mapped,bndryparam)/dx*len/2.0;

                        samplesv.at(g*2*halfpts + s) = point0; 
                        samplesv.at(g*2*halfpts + halfpts + s) = point1;


                    } else {
                        vertex point0 = mapped - (halfpts - 0.5 - s)*dx*unitNormal; 
                        vertex point1 = mapped + (0.5+s)*dx*unitNormal;

                        samples.at(g*2*halfpts + s) = diffBndry(mapped,bndryparam)/dx*len/2.0;
                        samples.at(g*2*halfpts + halfpts + s) = 
		  						 my_recon.at(FlatIndic(mi,cellid))->eval(locvals, point1, stenlg, stensm)/dx*len/2.0;

                        samplesv.at(g*2*halfpts + s) = point0;
                        samplesv.at(g*2*halfpts + halfpts + s) = point1;

	    	          }
                } else {
                    // Horizontal
                    if (cellid[1] == 0){
                        // Down sampling values should be replaced with dirichlet values
                        vertex point0 = mapped - (halfpts - 0.5 - s)*dx*unitNormal; 
                        vertex point1 = mapped + (0.5+(s-1))*dx*unitNormal;
    
                        samples.at(g*2*halfpts + s) = 
		    					 my_recon.at(FlatIndic(mi,cellid))->eval(locvals, point0, stenlg, stensm)/dx*len/2.0;
//                        samples.at(g*2*halfpts + halfpts + s) = my_recon.at(cellid)->eval(locvals, point1, stenlg, stensm);
                        samples.at(g*2*halfpts + halfpts + s) = diffBndry(mapped,bndryparam)/dx*len/2.0;

                        samplesv.at(g*2*halfpts + s) = point0; 
                        samplesv.at(g*2*halfpts + halfpts + s) = point1;

		              } else {
                        vertex point0 = mapped - (halfpts - 0.5 - s)*dx*unitNormal; 
                        vertex point1 = mapped + (0.5+s)*dx*unitNormal;

//                        samples.at(g*2*halfpts + s) = diffBndry(mapped,bndryparam)/dx*len/2.0;
//                        samples.at(g*2*halfpts + halfpts + s) = 
//								 my_recon.at(FlatIndic(mi,cellid))->eval(locvals, point1, stenlg, stensm)/dx*len/2.0;
                        samples.at(g*2*halfpts + halfpts + s) = diffBndry(mapped,bndryparam)/dx*len/2.0;
                        samples.at(g*2*halfpts + s) = 
								 my_recon.at(FlatIndic(mi,cellid))->eval(locvals, point1, stenlg, stensm)/dx*len/2.0;

                        samplesv.at(g*2*halfpts + s) = point0; 
                        samplesv.at(g*2*halfpts + halfpts + s) = point1;
	    	          }
	             }

				}
        } else {
            for (int s=0; s<halfpts; s++){
                    vertex point0 = mapped - (halfpts - 0.5 - s)*dx*unitNormal; 
                    vertex point1 = mapped + (0.5+s)*dx*unitNormal;

                    samples.at(g*2*halfpts + s) = 0.0;
                    samples.at(g*2*halfpts + halfpts + s) = 0.0;

                    samplesv.at(g*2*halfpts + s) = point0; 
                    samplesv.at(g*2*halfpts + halfpts + s) = point1;
            } 
		  }
	 }

    return 1;
}

int TransportVariable::EvaluateSamplesEdge(const MeshInfo& mi, 
                                           int xSize, int ySize, int offset, 
                                           int xMaxCell, int yMaxCell,
					                            double ** locvals){ 

    const valarray<double>& gwe = GaussWeightsEdge;

    // Initializing lagrangian interpolation 
    int degree  = gwe.size() + 1;
    int halfPts = std::ceil((degree+1)/2.0);
    int numPts  = halfPts * 2;

    // Define edge indicators
    int edgepos = 0;
    int edgeneg = 0;

    indice cellpos {0,0};
    indice cellneg {0,0};

    // Check if this edge is on the boundary or not
    bool onbndry = false;
    int bndrytype = 0;

    double len = 0.0;
    vertex unitNormal;

    for (int j=0; j<ySize; j++){
    for (int i=0; i<xSize; i++){

        getNeighbors(xMaxCell, yMaxCell, xSize, ySize, i, j, 
                     edgepos, edgeneg, cellpos, cellneg, onbndry);

        int dof = j*xSize + i;
        // Edge dof
        dof += offset;

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
            getUniformEdge(edge, uniformEdge, cellneg, xSize, ySize, xMaxCell, yMaxCell);
        } else {
            uniformEdge = edge;
        }

        len = length(uniformEdge);
        unitNormal = UnitNormal(uniformEdge, len);

//           cout <<  " The normal vector is " << unitNormal[0] <<  "  "  << unitNormal[1] << endl;

        // Extract surface area 
        double posh = sqrt(mi.cellArea.at(FlatIndic(mi, cellpos)));
        double negh = sqrt(mi.cellArea.at(FlatIndic(mi, cellneg)));

        double h = 2.0 * ((posh < negh) ? posh : negh);

        // Compute sample interval
        double dx = h /(double)(numPts - 1);

        // Extract sample points
        if (onbndry){
//            cout << "This is a boundary edge" << endl;
//            cout << "The cell id is "  << cellpos[0] << "  " << cellpos[1] << " and  " 
//					  << cellneg[0] << "  " << cellneg[1] << endl;
            EvaluateSamples(mi, cellneg, uniformEdge, unitNormal, dx, len, xSize, ySize, xMaxCell, yMaxCell, 
									 locvals, samplingp.at(dof), samplingv.at(dof));
			
        } else {
            EvaluateSamples(uniformEdge, unitNormal, dx, len, halfpts, 
								    FlatIndic(mi,cellneg), FlatIndic(mi,cellpos), locvals, samplingp.at(dof),
									 samplingv.at(dof));
        } 

	 }}

    return 1;
}
