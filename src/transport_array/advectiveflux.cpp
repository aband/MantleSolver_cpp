double advflux_edge(){


}

int advflux_all(const vector<reconstruction>& my_recon,
                double ** localvals,
					 vector<double>& advflux,
					 const MeshInfo& mi){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    advflux.clear();
    advflux.resize(M*(N+1) + N*(M+1));

    int pos = 0;
    int neg = 0;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    vector<double> quadwts;
    vector<vertex> quadpts;
    // Horizontal edge first M*(N+1)
    for (int j=0; j<N-1; j++){
        for (int i=0; i<M; i++){
            pos = j*M+i;
				neg = (j+1)*M+i;
            // Get edge quad points values
            quadwts.clear(); quadwts.resize(gwe.size());
            quadpts.clear(); quadpts.resize(gwe.size());

            advflux.at(neg) = ;    
        }
    }

    // Vertical edge second N*(M+1)

    return 1;
}
