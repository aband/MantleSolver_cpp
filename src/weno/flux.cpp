#include "../../include/flux.h"

solution LaxFriedrichsFlux(const WenoMesh*& wm, int pos, double t, WenoReconst**& wr, 
                           point target_point, point_index& target, 
                           double (*funcX)(valarray<double>& point, const vector<double>& param),
                           double (*funcY)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncY)(valarray<double>& point, const vector<double>& param)){

    solution flux;
    solution u_in, u_out;
    point_index  neighbor; 
    cell_corners corner;

    int totali = wm->M+2*wm->ghost;

    // Define Operation set
    int addi[4] = {-1,0,1,0};
    int addj[4] = {0,-1,0,1};

    int corner_rotate[4][2] = {{0,1},{0,0},{1,0},{1,1}};

    neighbor = {target[0]+addi[pos], target[1]+addj[pos]};
 
    corner.push_back(wm->lmesh[(target[1]+corner_rotate[pos][1])*totali + target[0]+corner_rotate[pos][0]]);
    corner.push_back(wm->lmesh[(target[1]+corner_rotate[(pos+1)%4][1])*totali + target[0]+corner_rotate[(pos+1)%4][0]]);

    double len = length(corner);
    point norm = UnitNormal(corner,len);

    // Transform target cell index to weno reconstruction index
    int index_in = (target[1]-wm->ghost+1)*(wm->M+2)+(target[0]-wm->ghost+1);
    int index_out = (neighbor[1]-wm->ghost+1)*(wm->M+2)+(neighbor[0]-wm->ghost+1);

    u_in = wr[index_in]->PointReconstruction(wm, target_point);

    u_out = wr[index_out]->PointReconstruction(wm, target_point);

//    cout << endl << "The edge position is : " << pos << endl;

//    cout << endl << "u_in is " << u_in << " u_out is " << u_out << endl;

    // Compute fluxand stabilization parameter
    flux = (funcX(target_point,{u_in}) + funcX(target_point,{u_out}))*norm[0] + (funcY(target_point,{u_in}) + funcY(target_point,{u_out}))*norm[1];

//    cout << endl << "The corresponding flux is : " << flux << endl;

    // Local Lax-Friedrichs scheme
    double alphaLF = max( fabs(dfuncX(target_point,{u_in})*norm[0] + dfuncY(target_point,{u_in})*norm[1]) , 
                          fabs(dfuncX(target_point,{u_in})*norm[0] + dfuncY(target_point,{u_in})*norm[1]) );

    // Global Lax-Friedrichs scheme
    alphaLF = 1.0;

    flux = 0.5 * (flux - alphaLF*(u_out-u_in));

    return flux;
}

solution TotalFlux(const WenoMesh*& wm, int pos, double t,
                   point_index& target, WenoReconst**& wr,
                   double (*funcX)(valarray<double>& point, const vector<double>& param),
                   double (*funcY)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                   double (*dfuncY)(valarray<double>& point, const vector<double>& param)){

    solution work = 0.0;
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    cell_corners corner;
    int totali = wm->M+2*wm->ghost;

    int corner_rotate[4][2] = {{0,1},{0,0},{1,0},{1,1}};
    corner.push_back(wm->lmesh[(target[1]+corner_rotate[pos][1])*totali + target[0]+corner_rotate[pos][0]]);
    corner.push_back(wm->lmesh[(target[1]+corner_rotate[(pos+1)%4][1])*totali + target[0]+corner_rotate[(pos+1)%4][0]]);

    double len = length(corner);

    for (int g=0; g<3; g++){
        valarray<double> mapped = GaussMapPointsEdge({gpe[g]},corner);
        work += gwe[g] * LaxFriedrichsFlux(wm, pos, t, wr, mapped, target, funcX, funcY, dfuncX, dfuncY);
    }
   
    work *= len/2.0; 

    return work;
}
