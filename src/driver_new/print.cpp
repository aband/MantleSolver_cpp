#include "driver.h"

static char * GetFilename(const char * fieldname, int mark){

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

int Driver::printGrid(){


    FILE * recongridx = fopen("recongridx.dat", "w");
    FILE * recongridy = fopen("recongridy.dat", "w");


    // Create sampling
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

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++) {

        for (int l=0; l<sampleSet.size(); l++){

            for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

               vertexSet corners = extractCorners(mi, {i,j});

               for (int g=0; g<sample1.size(); g++){

                   vertex mapped = GaussMapPointsFace(sampleSet.at(l)[g], corners); 

                   fprintf(recongridx, "%.16f ", mapped[0]);
                   fprintf(recongridy, "%.16f ", mapped[1]);
               }
            }
        }
    }

    fclose(recongridx);
    fclose(recongridy);

    return 1;
}

int Driver::printPhase(bool update, int mark){

    // Sometimes update reconstruction sometimes not  
	 // Print phase evaluation results at given data points
    FILE *fp = fopen(GetFilename("porosity", mark), "w");
    FILE *ft = fopen(GetFilename("temperature", mark), "w");
    FILE *fphase = fopen(GetFilename("phase", mark), "w");
    FILE *ftm= fopen(GetFilename("meltT", mark), "w");
    FILE *fopx = fopen(GetFilename("opx", mark), "w");
    FILE *fcd = fopen(GetFilename("cd", mark), "w");
    FILE *fhd = fopen(GetFilename("hd", mark), "w"); 

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    if (update){
        // Update reconstructions of CD
        vector<double> sigma_lg;
        sigma_lg.resize(stenlg.size());

        vector<double> sigma_sm;
        sigma_sm.resize(stensm.size());

        Vec localvecCD; 
        double ** locvalsCD;

        // Distribute local part to local vectors.
        PetscCall(DMGetLocalVector(dmu, &localvecCD)); 

        PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localvecCD));
        PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localvecCD));

        PetscCall(DMDAVecGetArray(dmu, localvecCD, &locvalsCD));

        for (int s=0; s<stenlg.size(); s++){
            sigma_lg.at(s) = stenlg.at(s).sigma(locvalsCD);
        }

        for (int s=0; s<stensm.size(); s++){
            sigma_sm.at(s) = stensm.at(s).sigma(locvalsCD);
        }

        // Setup nonlinear weights
        for (int s=0; s<my_recon_CD.size(); s++){
            my_recon_CD.at(s).extractsigma(sigma_lg, sigma_sm);
            my_recon_CD.at(s).setWgts(L_*H_/(double)M/(double)N);
        }

        // Update reconstructions of HD 
        sigma_lg.clear();
        sigma_lg.resize(stenlg.size());

        sigma_sm.clear();
        sigma_sm.resize(stensm.size());

        Vec localvecHD; 
        double ** locvalsHD;

        // Distribute local part to local vectors.
        PetscCall(DMGetLocalVector(dmu, &localvecHD)); 

        PetscCall(DMGlobalToLocalBegin(dmu, globalHD, INSERT_VALUES, localvecHD));
        PetscCall(DMGlobalToLocalEnd(dmu, globalHD, INSERT_VALUES, localvecHD));

        PetscCall(DMDAVecGetArray(dmu, localvecHD, &locvalsHD));

        for (int s=0; s<stenlg.size(); s++){
            sigma_lg.at(s) = stenlg.at(s).sigma(locvalsHD);
        }

        for (int s=0; s<stensm.size(); s++){
            sigma_sm.at(s) = stensm.at(s).sigma(locvalsHD);
        }

        // Setup nonlinear weights
        for (int s=0; s<my_recon_HD.size(); s++){
            my_recon_HD.at(s).extractsigma(sigma_lg, sigma_sm);
            my_recon_HD.at(s).setWgts(abs(L_*H_)/(double)M/(double)N);
        }

    }

    // Get local values HD
    Vec localHD;

    double ** lHD;

    PetscCall(DMGetLocalVector(dmu, &localHD));

    PetscCall(DMGlobalToLocalBegin(dmu, globalHD, INSERT_VALUES, localHD));
    PetscCall(DMGlobalToLocalEnd(dmu, globalHD, INSERT_VALUES, localHD));

    PetscCall(DMDAVecGetArray(dmu, localHD, &lHD));

    // Get local values CD

    Vec localCD;

    double ** lCD;

    PetscCall(DMGetLocalVector(dmu, &localCD));

    PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localCD));
    PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localCD));

    PetscCall(DMDAVecGetArray(dmu, localCD, &lCD));

    // Create sampling
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

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++) {

        for (int l=0; l<sampleSet.size(); l++){

            for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

               vertexSet corners = extractCorners(mi, {i,j});

               for (int g=0; g<sample1.size(); g++){

                   vertex mapped = GaussMapPointsFace(sampleSet.at(l)[g], corners); 

                   double lithoP = myPhase->pPtr->GetStaticP(mapped[1]*(-1), myPhase->pp->l0);

                   const int idx = FlatIndic(mi, {i,j});
                   double HD = my_recon_HD.at(idx).eval(lHD, mapped, stenlg, stensm);
                   double CD = my_recon_CD.at(idx).eval(lCD, mapped, stenlg, stensm);

                   fprintf(fhd, "%.16f ", HD);
                   fprintf(fcd, "%.16f ", CD);  

                   myPhase->pPtr->evalPhase(HD, CD, lithoP);

                   // Porosity
                   fprintf(fp, "%.16f ", myPhase->pPtr->pc.phil);

                   fprintf(fphase, "%d ", myPhase->pPtr->pc.region);

                   fprintf(ftm, "%.12f ", myPhase->pPtr->GetTDp(myPhase->pPtr->TDe0, lithoP));

                   fprintf(fopx, "%.12f ", myPhase->pPtr->pc.phi2);

               }
            }
        }
    }

    fclose(fp);
    fclose(ft);
    fclose(fphase);
    fclose(ftm);
    fclose(fopx);

    fclose(fhd);
	 fclose(fcd);

  
    return 1;
}
