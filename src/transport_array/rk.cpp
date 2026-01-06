#include "rk.h"

static int getflux(const MeshInfo& mi, Vec * innow, Vec * influx, DM dmu, DM dmmesh,
                   vector<reconstruction>& my_recon,
						 vector<tensorstencilpoly>& sten_lg,
						 vector<tensorstencilpoly>& sten_sm,
						 double t){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];
  
    Vec now  = *innow; 
    Vec flux = *influx;

    vector<double> sigma_lg;
    sigma_lg.resize(sten_lg.size());

    vector<double> sigma_sm;
    sigma_sm.resize(sten_sm.size());

    vector<double> advflux;

    Vec localu;

    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, now, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, now, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    double ** f;
    DMDAVecGetArray(dmu, flux, &f);

    for (int s=0; s<sten_lg.size(); s++){
        sigma_lg.at(s) = sten_lg.at(s).sigma(lu);
    }

    for (int s=0; s<sten_sm.size(); s++){
        sigma_sm.at(s) = sten_sm.at(s).sigma(lu);
    }

    // Setup nonlinear weights
    for (int s=0; s<my_recon.size(); s++){
        my_recon.at(s).extractsigma(sigma_lg, sigma_sm);
        my_recon.at(s).setWgts(1.0/(double)M/(double)N);
    }

    computeEdgeFlux(advflux, t, mi, lu, my_recon, sten_lg, sten_sm);  

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

       double area = mi.cellArea.at(FlatIndic(mi, {i,j}));

       int bottom = j*M+i;
       int top    = (j+1)*M+i;
       int left   = M*(N+1) + j*(M+1) + i;
       int right  = M*(N+1) + j*(M+1) + i+1;

       f[j][i] = (advflux.at(bottom) - advflux.at(top) + advflux.at(left) - advflux.at(right))/area;

    }}

    DMDAVecRestoreArray(dmu, flux, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}

int rk1(double dt, int Nt, Vec * insol, const MeshInfo& mi, 
        DM dmu, DM dmmesh,
        vector<reconstruction>& my_recon,
		  vector<tensorstencilpoly>& sten_lg,
		  vector<tensorstencilpoly>& sten_sm){

    Vec sol = *insol;
    int event = 1;

    for (int t=0; t<Nt; t++){

        Vec flux;
        VecDuplicate(sol, &flux);

        getflux(mi, &sol, &flux, dmu, dmmesh, my_recon, sten_lg, sten_sm, t*dt);

        VecAXPY(sol, -1*dt, flux);
    }

    return 1;
}

int rk2(double dt, int Nt, Vec * insol, const MeshInfo& mi, 
        DM dmu, DM dmmesh,
        vector<reconstruction>& my_recon,
		  vector<tensorstencilpoly>& sten_lg,
		  vector<tensorstencilpoly>& sten_sm){

    Vec sol = *insol;
    int event = 1;

    Vec temp;
    VecDuplicate(sol, &temp);
    VecCopy(sol, temp);

    for (int t=0; t<Nt; t++){

        Vec flux;
        VecDuplicate(sol, &flux);

        getflux(mi, &sol, &flux, dmu, dmmesh, my_recon, sten_lg, sten_sm, t*dt);

        VecAXPY(temp, -1*dt, flux);

        Vec flux2;
        VecDuplicate(sol, &flux2);

        getflux(mi, &temp, &flux2, dmu, dmmesh, my_recon, sten_lg, sten_sm, t*dt);

        VecScale(sol, 0.5);
        VecAXPY(sol, 0.5, temp);
        VecAXPY(sol, -0.5*dt, flux2);
    }

    return 1;
}
