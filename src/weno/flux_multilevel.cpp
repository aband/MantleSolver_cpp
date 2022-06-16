#include "../include/flux_multilevel.h"

double LaxFriedrichsFlux(const MeshInfo& mi, int pos, double t, WenoReconstruction**& wr,
                         point& target_point, point_index& target_index,
                         double (*funcX)(valarray<double>& point, const vector<double>& param),
                         double (*funcY)(valarray<double>& point, const vector<double>& param),
                         double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                         double (*dfuncY)(valarray<double>& point, const vector<double>& param)){

    double flux;
    double u_in, u_out;

    point_index neighbor;
    points_set corner;

    int fulllocalx = mi.localsize[0]+2*mi.ghost_vertx[0];

    // Define Operation set
    int addi[4] = {-1,0,1,0};
    int addj[4] = {0,-1,0,1};

    int corner_rotate[4][2] = {{0,1},{0,0},{1,0},{1,1}};

    neighbor = {target_index[0]+addi[pos], target_index[1]+addj[pos]};
 
    corner.push_back(mi.lmesh[(target_index[1]+corner_rotate[pos][1])*fulllocalx + 
                                target_index[0]+corner_rotate[pos][0]]);
    corner.push_back(mi.lmesh[(target_index[1]+corner_rotate[(pos+1)%4][1])*fulllocalx + 
                                target_index[0]+corner_rotate[(pos+1)%4][0]]);

    double len = length(corner);
    point norm = UnitNormal(corner,len);

    // Transform target cell index to weno reconstruction index
    int index_in = (target_index[1]-mi.ghost_vertx[1]+1)*(mi.localsize[0]+2)+
                   (target_index[0]-mi.ghost_vertx[0]+1);
    int index_out = (neighbor[1]-mi.ghost_vertx[1]+1)*(mi.localsize[0]+2)+
                    (neighbor[0]-mi.ghost_vertx[0]+1);

    u_in = wr[index_in]->PointValueReconstruction(mi, target_point);
    u_out = wr[index_out]->PointValueReconstruction(mi, target_point);

    // Compute fluxand stabilization parameter
    flux = (funcX(target_point,{u_in}) + funcX(target_point,{u_out}))*norm[0] + (funcY(target_point,{u_in}) + funcY(target_point,{u_out}))*norm[1];

//    cout << endl << "The corresponding flux is : " << flux << endl;

    // Local Lax-Friedrichs scheme
    //double alphaLF = max( fabs(dfuncX(target_point,{u_in})*norm[0] + dfuncY(target_point,{u_in})*norm[1]) , 
    //                     fabs(dfuncX(target_point,{u_in})*norm[0] + dfuncY(target_point,{u_in})*norm[1]) );

    // Global Lax-Friedrichs scheme
    double alphaLF = 1.0;

    flux = 0.5 * (flux - alphaLF*(u_out-u_in));

    return flux;
}

double TotalFlux(const MeshInfo& mi, int pos, double t,
                 point_index& target_index, WenoReconstruction**& wr,
                 double (*funcX)(valarray<double>& point, const vector<double>& param),
                 double (*funcY)(valarray<double>& point, const vector<double>& param),
                 double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                 double (*dfuncY)(valarray<double>& point, const vector<double>& param)){

    double work;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    points_set corner;

    int fulllocalx = mi.localsize[0]+2*mi.ghost_vertx[0];

    int corner_rotate[4][2] = {{0,1},{0,0},{1,0},{1,1}};
    corner.push_back(mi.lmesh[(target_index[1]+corner_rotate[pos][1])*fulllocalx + 
                                target_index[0]+corner_rotate[pos][0]]);
    corner.push_back(mi.lmesh[(target_index[1]+corner_rotate[(pos+1)%4][1])*fulllocalx + 
                                target_index[0]+corner_rotate[(pos+1)%4][0]]);

    double len = length(corner);

    for (int g=0; g<3; g++){
        valarray<double> mapped = GaussMapPointsEdge({gpe[g]},corner);
        work += gwe[g] * LaxFriedrichsFlux(mi,pos,t,wr,mapped,target_index,funcX,funcY,dfuncX,dfuncY); 
    }

    work *= len/2.0;

    return work;
}
