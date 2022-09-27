#include "../include/flux_multilevel.h"

double LaxFriedrichsFlux(const MeshInfo& mi, int pos, double t, vector<WenoReconstruction*>& wr,
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
                 point_index& target_index, vector<WenoReconstruction*>& wr,
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

// Create Derivative for a single Lax-Friedrich type flux for Jacobian
vector<double> DeriveLaxFriedrichFlux(const MeshInfo& mi, double t,
                                      point_index& target_index, vector<WenoReconstruction*>& wr,
                                      double (*funcX)(valarray<double>& point, const vector<double>& param),
                                      double (*funcY)(valarray<double>& point, const vector<double>& param),
                                      double (*dfuncX)(valarray<double>& point, const vector<double>& param),
                                      double (*dfuncY)(valarray<double>& point, const vector<double>& param)){
    
    vector<double> work;

    // Copy gauss points and weights
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    int stencilSizeX = wr[0]->GetStencilSizeX();
    int stencilSizeY = wr[0]->GetStencilSizeY();

    int derivativeSize = (stencilSizeX+2)*(stencilSizeY+2);

    work.resize(derivativeSize,0.0);

    int shiftSet[4][2] = {{-1,0},{0,-1},{1,0},{0,1}};

    int fulllocalx = mi.localsize[0]+2*mi.ghost_vertx[0];

    int corner_rotate[4][2] = {{0,1},{0,0},{1,0},{1,1}};

    // Define Operation set
    int addi[4] = {-1,0,1,0};
    int addj[4] = {0,-1,0,1};

    // Global Lax-Friedrichs scheme 
    double alphaLF = 1.0;

    // Loop through four edges
    for (int pos=0; pos<4; pos++){

        point_index  neighbor;
        neighbor = {target_index[0]+addi[pos], target_index[1]+addj[pos]};

        points_set corner;
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

        // Loop through gauss points
        for (int g=0; g<3; g++){
            valarray<double> mapped = GaussMapPointsEdge({gpe[g]},corner);
            vector<double> derivIn  = wr[index_in]->PseudoDerivativeWenoReconst(mi, mapped);
            vector<double> derivOut = wr[index_out]->PseudoDerivativeWenoReconst(mi, mapped);
            assert(derivIn.size() == derivOut.size());

            double u_in = wr[index_in]->PointValueReconstruction(mi, mapped);
            double u_out = wr[index_out]->PointValueReconstruction(mi, mapped);

            for (int s=0; s<derivIn.size(); s++){
                int localix = s%stencilSizeX;
                int localiy = s/stencilSizeX;

                int derivIndexIn  = (localiy + 1)*(stencilSizeX+2) + localix+1;
                int derivIndexOut = (localiy + 1 + shiftSet[pos][1])*(stencilSizeX+2) + 
                                     localix+1 + shiftSet[pos][0];

                work[derivIndexIn] += 0.5*derivIn[s]*(dfuncX(mapped,{u_in})*norm[0] + 
                                                      dfuncY(mapped,{u_in})*norm[1] +
                                                      alphaLF * 1.0)
                                      *gwe[g]*len/2.0;

                work[derivIndexOut] += 0.5*derivOut[s]*(dfuncX(mapped,{u_out})*norm[0] + 
                                                        dfuncY(mapped,{u_out})*norm[1] -
                                                        alphaLF * 1.0)
                                      *gwe[g]*len/2.0;

            }
        }

    }

    return work;
}
