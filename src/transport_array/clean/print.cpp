#include "transport.h"

int TransportVariable::Print(const MeshInfo& mi, const char * filename){

    FILE * sol = fopen(filename,"w");

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        for (int e=0; e<4; e++){
        for (int g=0; g<3; g++){

            fprintf(sol, "%.12f ", my_recon.at(j*mi.MPIglobalCellSize[0] + i)->elem_val.at(e*3+g));

        }} fprintf(sol, "\n");

    }fprintf(sol,"\n");}

    fclose(sol);

    return 1;
}

int TransportVariable::PrintSample(const MeshInfo& mi, int xSize, int ySize, int offset, 
                                   const char * filename, const char * filenamevx, const char * filenamevy){

    FILE * sample  = fopen(filename, "w");
    FILE * samplevx = fopen(filenamevx, "w");
    FILE * samplevy = fopen(filenamevy, "w");

    for (int j=0; j<ySize; j++){
    for (int i=0; i<xSize; i++){

        int dof = j*xSize + i + offset;

        for (int s=0; s<numpts; s++){
            fprintf(sample, "%.16f ", samplingp.at(dof).at(6+s));
            fprintf(samplevx, "%.16f ", samplingv.at(dof).at(6+s)[0]);
            fprintf(samplevy, "%.16f ", samplingv.at(dof).at(6+s)[1]);
        }

	 } fprintf(sample, "\n");
      fprintf(samplevx, "\n");
      fprintf(samplevy, "\n");}

    fclose(sample);
    fclose(samplevx);
    fclose(samplevy);

    return 1;
}

int TransportVariable::PrintEdgeSample(const MeshInfo& mi, int xSize, int ySize, int offset, 
                                       const char * filename, const char * filenamevx, const char * filenamevy){

    FILE * sample  = fopen(filename, "w");
    FILE * samplevx = fopen(filenamevx, "w");
    FILE * samplevy = fopen(filenamevy, "w");
    for (int j=0; j<ySize; j++){
    for (int i=0; i<xSize; i++){
        int dof = j*xSize + i + offset;
//		  int count = 0;
        for (int g=0; g<3; g++){
        for (int s=0; s<numpts; s++){
					 
            fprintf(sample, "%.16f ", samplingp.at(dof).at(g*numpts+s));
            fprintf(samplevx, "%.16f ", samplingv.at(dof).at(g*numpts+s)[0]);
            fprintf(samplevy, "%.16f ", samplingv.at(dof).at(g*numpts+s)[1]);
//				count ++;
        }}
	   fprintf(sample, "\n");
      fprintf(samplevx, "\n");
      fprintf(samplevy, "\n");
	 }}

    fclose(sample);
    fclose(samplevx);
    fclose(samplevy);

    return 1;
}

int TransportVariable::PrintPatch(const MeshInfo& mi, int m, int n, double ** locvals, 
                                  const char * filename, const char * filenamevx, const char * filenamevy){

    FILE * sol = fopen(filename, "w");

    FILE * vx = fopen(filenamevx, "w");
    FILE * vy = fopen(filenamevy, "w");

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    double dx = 2.0/(double) (m-1);
    double dy = 2.0/(double) (n-1);

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

        vertexSet corners = extractCorners(mi, {i,j});
        int dof = j*M+i;

        for (int jj=0; jj<n; jj++){
        for (int ii=0; ii<m; ii++){

            vertex sample {-1.0 + ii*dx, -1.0 + jj*dy};

            // Extract four corners
            vertex mapped = GaussMapPointsFace(sample, corners);
//cout <<"( " << mapped[0] << "  " << mapped[1] << ")   ";
            fprintf(sol, "%.16f ", my_recon.at(dof)->eval(locvals, mapped, stenlg, stensm));            
 
            fprintf(vx, "%.16f ", mapped[0]);            
            fprintf(vy, "%.16f ", mapped[1]);            
 
//        }cout << endl;}
        }}
        fprintf(sol, "\n");
        fprintf(vx, "\n");
        fprintf(vy, "\n");
//cout << endl;
	 }}

    fclose(sol);
    fclose(vx);
    fclose(vy);

    return 1;
}
