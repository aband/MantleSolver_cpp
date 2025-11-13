#include "integral.h"

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

valarray<double> GaussMapPointsFace(valarray<double> ref, 
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

valarray<double> GaussMapPointsEdge(valarray<double> ref,
                                    const vector< valarray<double> >& corner){
    assert(corner.size()==2);

    valarray<double> mid = (corner[0]+corner[1])/2.0;

    valarray<double> temp = (corner[0]-corner[1]);
    temp = temp * temp;
    double len = sqrt(temp.sum());

    valarray<double> mapped = mid + abs(temp)/len*ref[0]/2;

    return mapped;
}

double GaussJacobian(valarray<double> ref,
                     const vector< valarray<double> >& corner){
    assert(corner.size()==4);

    double jac = 1.0/16 * ( (- corner[0][0]*(1.0-ref[1]) + corner[1][0]*(1.0-ref[1]) + corner[2][0]*(1.0+ref[1]) - corner[3][0]*(1.0+ref[1]) ) *
                            (- corner[0][1]*(1.0-ref[0]) - corner[1][1]*(1.0+ref[0]) + corner[2][1]*(1.0+ref[0]) + corner[3][1]*(1.0-ref[0]) ) -
                            (- corner[0][0]*(1.0-ref[0]) - corner[1][0]*(1.0+ref[0]) + corner[2][0]*(1.0+ref[0]) + corner[3][0]*(1.0-ref[0]) ) *
                            (- corner[0][1]*(1.0-ref[1]) + corner[1][1]*(1.0-ref[1]) + corner[2][1]*(1.0+ref[1]) - corner[3][1]*(1.0+ref[1]) ) );

    return jac; 
}

valarray<double> UnitNormal(const vector< valarray<double> >& corner,double len){
    /*
     * The function calculates the unit outer normal vector on a given edge/face(later)
     * Rotation 90 degrees clockwise
     * Resulting in four outward unit normal vectors
     */
    valarray<double> work;
   
    work = corner[1]-corner[0];
    work = work.cshift(1);
    work[1] = -work[1];
    work = work/len;

    return work;
}

double length(const vector< valarray<double> >& corner){

    valarray<double> vec = corner[1]-corner[0];
    vec *= vec;
    return sqrt(vec.sum()); 
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
                       double (*func)(const valarray<double>& point, const vector<double>& param)){

    assert(corner.size() == 4);

    double work = 0.0;

    vector< valarray<double> > tmp = corner;
    /*
     *Transform original corner coordinates with given
     *parameter h and center point. If no transform, pass
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

double NumIntegralFace(const vector< valarray<double> >& corner, const vector<int>& param,
                       const valarray<double>& center, const double& h, 
                       double (*func)(const valarray<double>& point, const vector<int>& param)){

    assert(corner.size() == 4);
 
    double work = 0.0;

    vector< valarray<double> > tmp = corner;
    /*
     *Transform original corner coordinates with given
     *parameter h and center point. If no transform, pass
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

double NumIntegralEdge(const vector< valarray<double> >& corner, const vector<int>& param,
                       double (*funcX)(const valarray<double>& point, const vector<int>& param), 
                       double (*funcY)(const valarray<double>& point, const vector<int>& param)){

    assert(corner.size() == 2);

    double work = 0.0;
     
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(corner);

    valarray<double> norm = UnitNormal(corner,len);

    for (int g=0; g< 3; g++){
        valarray<double> mapped = GaussMapPointsEdge({gpe[g]},corner);
        printf("The mapped gauss points on the edge are (%f,%f) \n",mapped[0],mapped[1]);
        work += gwe[g] * (funcX(mapped,param)*norm[0] + funcY(mapped,param)*norm[1]);
    }
 
    work *= len/2.0; 

    return work;
}

double NumIntegralEdge(const vector< valarray<double> >& corner, const vector<double>& param,
                       double (*funcX)(const valarray<double>& point, const vector<double>& param), 
                       double (*funcY)(const valarray<double>& point, const vector<double>& param)){

    assert(corner.size() == 2);

    double work = 0.0;
     
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    double len = length(corner);

    valarray<double> norm = UnitNormal(corner,len);

    for (int g=0; g< 3; g++){
        valarray<double> mapped = GaussMapPointsEdge({gpe[g]},corner);
        work += gwe[g] * (funcX(mapped,param)*norm[0] + funcY(mapped,param)*norm[1]);
    }
 
    work *= len/2.0; 

    return work;
}
