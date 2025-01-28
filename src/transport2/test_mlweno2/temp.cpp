#include "temp.h"

double func(const vertex& point,
            const vector<double>& param){

    // Initial condition

    return point[0]*point[0];
}

int RK(double dt, int Nt, Vec * insol, const MeshInfo& mi, multilevel& ml, const mluse& use, DM dmu, DM dmmesh){

    Vec sol  = *insol;
    Vec flux;
    VecDuplicate(sol, &flux);

    for (int t=0 ; t<Nt; t++){

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




    // Loop through physical domain
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        f[j][i] += vertedgeflux({i,j}) - vertedgeflux({i+1,j}) + horiedgeflux({i,j}) - horiedgeflux({i,j+1});

    }}

    DMDAVecRestoreArray(dmu, flux, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}

int advfunc(const vector<double>& uin, const vector<double>& uout,
                  vector<double>& fin,       vector<double>& fout,
            const vector<vertex>& param,     vector<double>& LF,
            const vertex& unitnormal){

    // Representing transport with prescribed velocity
    // u_t + v u = 0;

    for (int i=0; i<uin.size(); i++){
        double vel = param.at(i)[0] * unitnormal[0] + param.at(i)[1] * unitnormal[1];
        fin[i]  = vel*uin.at(i); 
        fout[i] = vel*uout.at(i); 
        LF[i]   = vel;
    }

    return 1;
}


