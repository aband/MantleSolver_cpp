#ifndef FUNC_H_
#define FUNC_H_

#include <valarray>
#include <vector>

using namespace std;

double funcX(valarray<double>& target, const vector<double>& param);
double funcY(valarray<double>& target, const vector<double>& param);
double dfuncX(valarray<double>& target, const vector<double>& param);
double dfuncY(valarray<double>& target, const vector<double>& param);

double Initial_Condition(valarray<double>& target, const vector<double>& param);

#endif
