#ifndef FUNC_H_
#define FUNC_H_

#include <valarray>
#include <vector>

using namespace std;

double funcF(double t, valarray<double>& target, double u);
double funcG(double t, valarray<double>& target, double u);
double dfuncF(double t, valarray<double>& target, double u);
double dfuncG(double t, valarray<double>& target, double u);

double Initial_Condition(valarray<double>& target, const vector<double>& param);

#endif
