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

/*
 *Define numerical integral function with function
 *overloading in c++.
 */
double NumIntegralFace(const vector< valarray<double> >& corner, const vector<double>& param,
                       const valarray<double>& center, const double& h, 
                       double (*func)(valarray<double>& point, const vector<double>& param));

double NumIntegralFace(const vector< valarray<double> >& corner, const vector<int>& param,
                       const valarray<double>& center, const double& h, 
                       double (*func)(valarray<double>& point, const vector<int>& param));

double NumIntegralEdge(const vector< valarray<double> >& corner, const vector<int>& param,
                       double (*funcX)(valarray<double>& point, const vector<int>& param),
                       double (*funcY)(valarray<double>& point, const vector<int>& param));

#endif
