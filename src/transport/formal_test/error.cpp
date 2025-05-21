// Returns reconstruction error computed in selected norm.
#include "error.h"

std::string location(const MeshInfo& mi, const indice& gcell){

    //if (gcell[1] == 0 || gcell[1] == mi.MPIglobalCellSize[1]-1){
    //    return "side";
    //} else {
    //    return "interior";
   // }

    return "interior";
}

int reconPlot(const MeshInfo& mi, multilevel& ml, mluse& use, int mark, Vec * global, bool grid, DM dmu, double h0){

    Vec temp = *global;

    const char * fieldname = "reconSol";

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];

    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename,"w");

    FILE * gridreconx = fopen("gridreconx.dat", "w"); 
    FILE * gridrecony = fopen("gridrecony.dat", "w");

    vector<vertex> sample1 = {{-1+1e-3,-1+1e-3},
                              { 0     ,-1+1e-3},
                              { 1-1e-3,-1+1e-3}};

    vector<vertex> sample2 = {{-1+1e-3, 0},
                              { 0     , 0},
                              { 1-1e-3, 0}};

    vector<vertex> sample3 = {{-1+1e-3, 1-1e-3},
                              { 0     , 1-1e-3},
                              { 1-1e-3, 1-1e-3}};

    vector<vector<vertex>> sampleSet = {sample1,sample2,sample3};

    Vec localvec;
    double ** locvals;

    PetscCall(DMGetLocalVector(dmu, &localvec)); 

    PetscCall(DMGlobalToLocalBegin(dmu, temp, INSERT_VALUES, localvec));
    PetscCall(DMGlobalToLocalEnd(dmu, temp, INSERT_VALUES, localvec));

    PetscCall(DMDAVecGetArray(dmu, localvec, &locvals));

    ml.updatesigma(locvals);
    Tensor<weights> allwgts;
 
    use.computeWgts(ml, mi, h0, allwgts, location);
 
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        // loop through scanning levels
        for (int l=0; l<3; l++){

            for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

                // Extract corners of the selected cell
                vertexSet corners = extractCorners(mi, {i,j});

                for (int g=0; g<3; g++){
                    vertex mapped = GaussMapPointsFace(sampleSet.at(l)[g], corners);

                    double val = use.eval(mapped, ml, location(mi, {i,j}), allwgts({i,j}), {i,j}, locvals); 

                    fprintf(sol, "%.12f ", val);

                    if (grid){
                        fprintf(gridreconx, "%.12f ",mapped[0]); 
                        fprintf(gridrecony, "%.12f ",mapped[1]);
                    }
                }

            }

        }
    }

    DMDAVecRestoreArray(dmu,localvec,&locvals);
    DMRestoreLocalVector(dmu, &localvec); 

    fclose(gridreconx);
    fclose(gridrecony);

    fclose(sol);

    return 1;
}

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

    use.computeWgts(ml, mi, h0, allwgts,location);

    // Update edgeflux
    Tensor<double> horiedgeflux = Tensor<double>(2);
    horiedgeflux.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgeflux = Tensor<double>(2);
    vertedgeflux.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    updateEdgeFlux(vertedgeflux, horiedgeflux, mi, lu, use, ml, allwgts);

    // Loop through physical domain
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        f[j][i] = getcellflux(mi, {i,j}, vertedgeflux, horiedgeflux);
    }}

    DMDAVecRestoreArray(dmu, flux, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}

int simpleRK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml,  mluse& use, DM dmu, DM dmmesh){

    Vec sol  = *insol;
    int event = 1;
    printSol(event,&sol,mi);
    for (int t=0 ; t<Nt; t++){

        Vec flux;
        VecDuplicate(sol, &flux);

        getflux(mi, ml, use, &sol, &flux, dmu, dmmesh);

        VecAXPY(sol, -1*dt, flux);

        if (t%5 == 0){
        event ++;
        printSol(event,&sol,mi);
        }
    }
    printSol(event,&sol,mi);

    return 1;
}

int SSP3RK(){


    return 1;
}
