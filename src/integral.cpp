#include "../include/integral.h"

/*
 *Gauss-Legendre points and corresponding weights are defined as
 *external variables. 
 */
valarray<double> GaussWeightsEdge {5.0/9.0,8.0/9.0,5.0/9.0};

valarray<double> GaussPointsEdge {-sqrt(3.0/5.0),0.0,sqrt(3.0/5.0)};

valarray<double> GaussWeightsFace {25.0/81.0,40.0/81.0,25/81.0,
                                  40.0/81.0,64.0/81.0,40/81.0,
                                  25.0/81.0,40.0/81.0,25/81.0};

vector< valarray<double> > GaussPointsFace  { {-sqrt(3.0/5.0), -sqrt(3.0/5.0)}, 
                                              {0.0           , -sqrt(3.0/5.0)}, 
                                              {sqrt(3.0/5.0) , -sqrt(3.0/5.0)}, 
                                              {-sqrt(3.0/5.0), 0.0}, 
                                              {0.0           , 0.0}, 
                                              {sqrt(3.0/5.0) , 0.0},
                                              {-sqrt(3.0/5.0), sqrt(3.0/5.0)},
                                              {0.0           , sqrt(3.0/5.0)},
                                              {sqrt(3.0/5.0) , sqrt(3.0/5.0)} };

inline valarray<double> GaussMapPointsFace(valarray<double> ref, 
                                           const vector< valarray<double> >& corner){
    assert(corner.size()==4);

    valarray<double> mapped = {0.0,0.0};
    double N[4] = { 0.25*(1.0-ref[0])*(1-ref[1]), 0.25*(1.0+ref[0])*(1-ref[1]),
                    0.25*(1.0+ref[0])*(1+ref[1]), 0.25*(1.0-ref[0])*(1+ref[1]) };

    for (size_t i=0; i<corner.size(); i++){
        mapped += corner[i]*N[i];
    }

    return mapped;
}

inline valarray<double> GaussMapPointsEdge(valarray<double> ref,
                                           const vector< valarray<double> >& corner){
    assert(corner.size()==2);

    valarray<double> mid = (corner[0]+corner[1])/2.0;

    valarray<double> temp = (corner[0]-corner[1]);
    temp = temp * temp;
    double len = sqrt(temp.sum());

    valarray<double> mapped = mid + abs(temp)/len*ref[0];
     
    return mapped;
}

inline double GaussJacobian(valarray<double> ref,
                            const vector< valarray<double> >& corner){
    assert(corner.size()==4);

    double jac = 1.0/16 * ( (- corner[0][0]*(1.0-ref[1]) + corner[1][0]*(1.0-ref[1]) + corner[2][0]*(1.0+ref[1]) - corner[3][0]*(1.0+ref[1]) ) *
                            (- corner[0][1]*(1.0-ref[0]) - corner[1][1]*(1.0+ref[0]) + corner[2][1]*(1.0+ref[0]) + corner[3][1]*(1.0-ref[0]) ) -
                            (- corner[0][0]*(1.0-ref[0]) - corner[1][0]*(1.0+ref[0]) + corner[2][0]*(1.0+ref[0]) + corner[3][0]*(1.0-ref[0]) ) *
                            (- corner[0][1]*(1.0-ref[1]) + corner[1][1]*(1.0-ref[1]) + corner[2][1]*(1.0+ref[1]) - corner[3][1]*(1.0+ref[1]) ) );

    return jac; 
}

/*
 *The numbering order for each cell/element is fixed as
 *  3_______2
 *  |       |
 *  |       |
 *  |       |
 *  0_______1
 */

double NumIntegralFace(const vector< valarray<double> >& corner, const vector<double>& param,
                       const valarray<double>& center, const double& h, 
                       double (*func)(valarray<double>& point, const vector<double>& param)){
    double work = 0.0;

    vector< valarray<double> > tmp = corner;
    /*
     *Transform original corner coordinates with given
     *parameter h and center point. If transform, pass
     *in h=1.0 and center point as (0.0,0.0).
     */
    for (auto & p : tmp){
        p -= center;
        p = p/h; 
    }

    // Gauss Quadrature
    const valarray<double>& gwf = GaussWeightsFace;
    const vector< valarray<double> >& gpf = GaussPointsFace;

    for (size_t i=0; i<gpf.size(); i++){
        valarray<double> mapped = GaussMapPointsFace(gpf[i],tmp);
        double jac = abs(GaussJacobian(gpf[i],tmp));
        double gw = gwf[i];
        work += jac*gw*(*func)(mapped,param);
    }

    return work;
}

/*
 *Edge Integral will be added later
 */

/*
 *double NumIntegralEdge(const valarray<double>){
 *
 *
 *}
 */
