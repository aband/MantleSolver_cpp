#include "error.h"

int printexactsol(const MeshInfo& mi, double t, 
                  double (*func)(const vertex& point,
                                 const vector<double>& param), 
				      int mark, bool grid, const vector<double>& param){

    const char * fieldname = "exactsol";

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];

    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename,"w");

    FILE * exactgridx = fopen("exactgridx.dat", "w");
    FILE * exactgridy = fopen("exactgridy.dat", "w");

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

    // Print exact solution on the given sample points
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        // loop through scanning levels
        for (int l=0; l<3; l++){

            for (int i=0; i<mi.MPIglobalCellSize[0]; i++){
                vertexSet corners = extractCorners(mi, {i,j});

                for (int g=0; g<3; g++){
                    vertex mapped = GaussMapPointsFace(sampleSet.at(l)[g], corners);
                    fprintf(sol, "%.12f ", func(mapped, {t}));

                    if (grid) {
                        fprintf(exactgridx, "%.12f ",mapped[0]); 
                        fprintf(exactgridy, "%.12f ",mapped[1]);
                    }
                }
            }
        }
    } 

    return 1;
}


