#include "temp.h"

int printSol(int mark, Vec * global, const MeshInfo& mi){

    Vec temp = *global;

    const char * fieldname = "sol";

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename,"w");

    for(int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for(int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double val;
        indice index {i,j};

        int nelem = FlatIndic(mi, index);

        PetscCall(VecGetValues(temp,1, &nelem, &val));

        fprintf(sol, "%e ", val);

    }fprintf(sol, "\n");}

    fclose(sol);

    return 1;
}

int printGrid(const MeshInfo& mi){

    FILE *gridPorox = fopen("gridCellX.dat", "w");
    FILE *gridPoroy = fopen("gridCellY.dat", "w");

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        vertex local {0.0,0.0};

        vector<vertex> corners = extractCorners(mi, {i,j});

        vertex global = GaussMapPointsFace(local, corners);

        fprintf(gridPorox,"%f ",global[0]);
        fprintf(gridPoroy,"%f ",global[1]);

    }
    fprintf(gridPorox, "\n");
    fprintf(gridPoroy, "\n");}

    fclose(gridPorox);
    fclose(gridPoroy);

    return 1;
}

// Initial consition
double func(const vertex& point,
            const vector<double>& param){

    // Initial condition

    // Initialize with simple Reimann shock and rarefaction function
    // time inputed as param[0] 

    if (point[0] < 0.5 || point[0] >=(0.5*param[0]+1.5)){
        return 0;
    } else if (point[0]>=0.5 && point[0]<param[0]+0.5){
        return (point[0]-0.5)/param[0];
    } else if (point[0]>=point[0]+0.5 || point[0] <0.5*param[0]+1.5){
        return 1;
    } else {
        return 0;
    }

}

// Burgers for testing
int advfunc(const vector<double>& uin, const vector<double>& uout,
                  vector<double>& fin,       vector<double>& fout,
            const vector<vertex>& param,     vector<double>& LF,
            const vertex& unitnormal){

/*
    // Representing transport with prescribed velocity
    // u_t + v u = 0;
    for (int i=0; i<uin.size(); i++){
        double vel = param.at(i)[0] * unitnormal[0] + param.at(i)[1] * unitnormal[1];
        fin[i]  = vel*uin.at(i); 
        fout[i] = vel*uout.at(i); 
        LF[i]   = vel;
    }
*/

    for (int i=0; i<uin.size(); i++){
        fin[i]  = uin.at(i)*uin.at(i)/2.0; 
        fout[i] = uout.at(i)*uout.at(i)/2.0; 
        LF[i]   = 1.0;
    }

    return 1;
}

int RK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml,  mluse& use, DM dmu, DM dmmesh){

    Vec sol  = *insol;
    for (int t=0 ; t<Nt; t++){

        Vec flux;
        VecDuplicate(sol, &flux);

        getflux(mi, ml, use, &sol, &flux, dmu, dmmesh);

        VecAXPY(sol, dt, flux);

        if (t%5 == 0){
        printSol(t/5+1,&sol,mi);
        }
    }

    return 1;
}

int getflux(const MeshInfo& mi, multilevel& ml, mluse& use, Vec * innow, Vec * influx, DM dmu, DM dmmesh){

    Vec now  = *innow; 
    Vec flux = *influx;

    Vec localu;

    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, now, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, now, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    double ** f;
    DMDAVecGetArray(dmu, flux, &f);

    // Update non linear weights with current cell-averaged solution
    ml.updatesigma(lu);

    Tensor<weights> allwgts;
    double h0 = sqrt((mi.L*mi.H)/(double)(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]));
    use.computeWgts(ml, mi, h0, allwgts);

    // Update edgeflux
    Tensor<double> horiedgeflux = Tensor<double>(2);
    horiedgeflux.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgeflux = Tensor<double>(2);
    vertedgeflux.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    getedgefluxall(mi, ml, use, lu, allwgts, horiedgeflux, vertedgeflux);

    // Loop through physical domain
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        f[j][i] += vertedgeflux({i,j}) - vertedgeflux({i+1,j}) + horiedgeflux({i,j}) - horiedgeflux({i,j+1});

    cout << setw(6) << "At cell (" << i << ", " << j << ")" << endl;
    cout << setw(6) << std::right << std::scientific
         << "Left edge flux   : " << vertedgeflux({i,j})   << "  "
         << "Right edge flux  : " << vertedgeflux({i+1,j})  << "  "
         << "Bottom edge flux : " << horiedgeflux({i,j}) << "  "
         << "Top edge flux    : " << horiedgeflux({i,j+1})    << endl << endl;

    //    cout << f[j][i] << " " ; 
    }cout << endl;}

    DMDAVecRestoreArray(dmu, flux, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}

inline bool onboundary(int i, int j, const MeshInfo& mi, int flag){

    if (flag == 1){
        // Vertical
        if (i==0 || i==mi.MPIglobalCellSize[0]+1){
            return true;
        } else {
            return false;
        }


    } else {
        // horizontal
        if (j==0 || j==mi.MPIglobalCellSize[1]+1){
            return true;
        } else {
            return false;
        }

    }

}

int getedgefluxall(const MeshInfo& mi, multilevel& ml, mluse& use, double ** lu,
                   const Tensor<weights>& allwgts,
                   Tensor<double>& horiedgeflux, Tensor<double>& vertedgeflux){

    // Full serial, not intended for parallel
    // Used for implicit testing

    vertex unitnormal;
    double len;
    vector<double> uin;
    vector<double> uout;
    vector<vertex> param {{1.0,0.0}};

    indice gCellOut, gCellIn;
    edgeEnds<vertex> edgeEndsVertex;
    edgeEnds<indice> edgeEndsIndice;

    const valarray<double>& gpe = GaussPointsEdge;
    uin.resize(gpe.size());
    uout.resize(gpe.size());

    // Loop through vertical edges
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
        for (int i=0; i<mi.MPIglobalCellSize[0]+1; i++){

            double flux = 0.0;

            if (i==0 || i == mi.MPIglobalCellSize[0]){
                flux = fluxintegralbndry();
            } else {

                indice edgeindex {i,j};
                extractVertEdgeInfo(mi, edgeindex, mi.ghostShiftVertex, gCellOut, gCellIn, edgeEndsVertex, edgeEndsIndice);

                array<vertex, 2> edge {edgeEndsVertex.start, edgeEndsVertex.end};

                len = getEdgeLength(edge);

                unitnormal = getUnitNormal(edge, len);

                for (int g=0; g<gpe.size(); g++){
                    vertex mapped = GaussMapPointsEdge({gpe[g]}, {edge[0],edge[1]});
                    uin.at(g)  = use.eval(mapped, ml, "all", allwgts({gCellIn[0] , gCellIn[1] }), gCellIn,  lu);
                    uout.at(g) = use.eval(mapped, ml, "all", allwgts({gCellOut[0], gCellOut[1]}), gCellOut, lu);
                } 

                flux = fluxintegral(unitnormal, len, uin, uout, param); 
            }

            vertedgeflux({i,j}) = flux;        
        }
    }

    // Loop through horizontal edges
    for (int j=0; j<mi.MPIglobalCellSize[1]+1; j++){
        for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

            double flux = 0.0;

            if (j==0 || j==mi.MPIglobalCellSize[1]){
                flux = fluxintegralbndry();
            } else {

                indice edgeindex {i,j};
                extractHoriEdgeInfo(mi, edgeindex, mi.ghostShiftVertex, gCellOut, gCellIn, edgeEndsVertex, edgeEndsIndice);

                array<vertex,2> edge {edgeEndsVertex.start, edgeEndsVertex.end};

                len = getEdgeLength(edge);

                unitnormal = getUnitNormal(edge, len);

                for (int g=0; g<gpe.size(); g++){
                    vertex mapped = GaussMapPointsEdge({gpe[g]}, {edge[0],edge[1]});
                    uin.at(g)  = use.eval(mapped, ml, "all", allwgts({gCellIn[0] , gCellIn[1] }), gCellIn,  lu);
                    uout.at(g) = use.eval(mapped, ml, "all", allwgts({gCellOut[0], gCellOut[1]}), gCellOut, lu);
                } 

                flux = fluxintegral(unitnormal, len, uin, uout, param); 
            }

            horiedgeflux({i,j}) = flux;        

        }
    }

    return 1;
}
