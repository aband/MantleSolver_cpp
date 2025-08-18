#include "error.h"

double LnormError(){

    double work = 0.0;


    return work;
}

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
//    vector<vector<vertex>> sampleSet = {sample1,sample3};

    // Print exact solution on the given sample points
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        // loop through scanning levels
        for (int l=0; l<sampleSet.size(); l++){

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

int printreconsol(const vector<reconstruction>& my_recon, int M, int N, int mark){

    const char * fieldname = "reconSol";

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];

    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename,"w");

    for (int j=0; j<N; j++){
        for (int l=0; l<3; l++){
            for (int i=0; i<M; i++){
                for (int g=0; g<3; g++){
                    fprintf(sol, "%.12f ", my_recon.at(j*M+i).elem_val.at(l*3+g) );
                }
            }
        }
    }

    fclose(sol);

    return 1;
}

int printreconsol(vector<reconstruction>& my_recon, int M, int N, int mark, 
                  const vector<tensorstencilpoly>& sten_lg,
						const vector<tensorstencilpoly>& sten_sm,
						const MeshInfo& mi,
						double ** locvals){

    vector<vertex> sample = {{-1+1e-3,-1+1e-3},
                             { 0     ,-1+1e-3},
                             { 1-1e-3,-1+1e-3},
									  {-1+1e-3, 0},
							        { 0     , 0},
							        { 1-1e-3, 0},
									  {-1+1e-3, 1-1e-3},
                             { 0     , 1-1e-3},
                             { 1-1e-3, 1-1e-3}};

    vector<vertex> mapped; 
    mapped.resize(sample.size());

    // Evalutation at sample points
    for (int j=0; j<N; j++){

        for (int i=0; i<M; i++){

            // Extract corners of the selected cell
            vertexSet corners = extractCorners(mi, {i,j});

            for (int g=0; g<sample.size(); g++){
                mapped.at(g) = GaussMapPointsFace(sample.at(g), corners);
            }

            my_recon.at(j*M+i).eval(locvals, mapped, sten_lg, sten_sm);
            //if (j==midN && i==midM){
            //    my_recon.at(j*M+i).eval(locvals, mapped, sten5, sten3);
            //}
        }
    }

    return printreconsol(my_recon, M ,N, 1);
}

int printreconsol2(vector<reconstruction>& my_recon, int M, int N, int mark, 
                  const vector<tensorstencilpoly>& sten_lg,
						const vector<tensorstencilpoly>& sten_sm,
						const MeshInfo& mi,
						double ** locvals){

    const char * fieldname = "reconSol";

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];

    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename,"w");

//    FILE * exactgridx = fopen("exactgridx.dat", "w");
//    FILE * exactgridy = fopen("exactgridy.dat", "w");

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

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        // loop through scanning levels
        for (int l=0; l<sampleSet.size(); l++){

            for (int i=0; i<mi.MPIglobalCellSize[0]; i++){
                vertexSet corners = extractCorners(mi, {i,j});

                for (int g=0; g<3; g++){
                    vertex mapped = GaussMapPointsFace(sampleSet.at(l)[g], corners);
                    fprintf(sol, "%.12f ", my_recon.at(j*M+i).eval(locvals, mapped, sten_lg, sten_sm) );
                }
            }
        }
    } 

    return 1;

}
