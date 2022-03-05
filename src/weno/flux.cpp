#include "../../include/flux.h"

solution LaxFriedrichsFlux(const WenoMesh*& wm, solution u_in, int pos, double t,
                           index_set& StencilLarge, vector<index_set>& StencilSmall, point_index& target, 
                           vector<int>& Sorder, vector<int>& Lorder, point target_point,
                           double (*funcX)(valarray<double>& point, const vector<double>& param),
                           double (*funcY)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                           double (*dfuncY)(valarray<double>& point, const vector<double>& param)){

    solution flux;
    solution u_out;
    point_index neighbor;
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

    u_out = WenoPointReconst(StencilLarge,StencilSmall,wm,neighbor,Sorder,Lorder,target_point); 

    flux = (funcX(target_point,{u_in}) + funcX(target_point,{u_out}))*norm[0] + (funcY(target_point,{u_in}) + funcY(target_point,{u_out}))*norm[1];

    // Determine stabilization parameter
    double alphaLF = max( pow( pow(dfuncX(target_point,{u_in}),2)  + pow(dfuncY(target_point,{u_in}),2 ),0.5) , 
                          pow( pow(dfuncX(target_point,{u_out}),2) + pow(dfuncY(target_point,{u_out}),2),0.5) );

    flux = 0.5 * (flux - alphaLF*(u_out-u_in));

    return flux;
}
