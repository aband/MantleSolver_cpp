#ifndef INTEGRAL_H_
#define INTEGRAL_H_

#include <vector>
#include <valarray>
#include <iostream>
#include <cmath>
#include <assert.h>

using namespace std;

/*
 *Define external variables
 */
extern valarray<double> GaussWeightsEdge;

extern valarray<double> GaussPointsEdge;

extern valarray<double> GaussWeightsFace;

extern vector< valarray<double> > GaussPointsFace;

valarray<double> UnitNormal(const vector< valarray<double> >& corner, double len);

double length(const vector< valarray<double> >& corner);

valarray<double> GaussMapPointsEdge(valarray<double> ref, const vector< valarray<double> >& corner);

double GaussJacobian(valarray<double> ref, const vector< valarray<double> >& corner);
 
valarray<double> GaussMapPointsFace(valarray<double> ref, const vector< valarray<double> >& corner);
 
/*
 *Define numerical integral function with function
 *overloading in c++.
 */
double NumIntegralFace(const vector< valarray<double> >& corner, const vector<double>& param,
                       const valarray<double>& center, const double& h, 
                       double (*func)(const valarray<double>& point, const vector<double>& param));

double NumIntegralFace(const vector< valarray<double> >& corner, const vector<int>& param,
                       const valarray<double>& center, const double& h, 
                       double (*func)(const valarray<double>& point, const vector<int>& param));

double NumIntegralEdge(const vector< valarray<double> >& corner, const vector<int>& param,
                       double (*funcX)(const valarray<double>& point, const vector<int>& param),
                       double (*funcY)(const valarray<double>& point, const vector<int>& param));

double NumIntegralEdge(const vector< valarray<double> >& corner, const vector<double>& param,
                       double (*funcX)(const valarray<double>& point, const vector<double>& param),
                       double (*funcY)(const valarray<double>& point, const vector<int>& param));
#endif
