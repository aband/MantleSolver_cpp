#include "temp.h"

double func(const vertex& point,
            const vector<double>& param){

    // Initial condition

    return point[0]*point[0];
}

int RK(double dt, int Nt, Vec * insol, const MeshInfo& mi, const multilevel& ml, const mluse& use){

    Vec sol  = *insol;
    Vec flux;
    VecDuplicate(sol, &flux);

    for (int t=0 ; t<Nt; t++){

    }

    return 1;
}

int getflux(const MeshInfo& mi, const multilevel& ml, const mluse& use, Vec * innow, Vec * influx){

    Vec now  = *innow; 
    Vec flux = *influx;

    Vec localu;

    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, U, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, U, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    double ** f;
    DMDAVecGetArray(mi.dmu, F, &f);

    use.updatesigma(ml,lu);

    DMDAVecRestoreArray(mi.dmu, F, &f);
    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);


    return 1;
}
