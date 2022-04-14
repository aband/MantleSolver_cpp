#ifndef FUNC_H_
#define FUNC_H_

#include <valarray>
#include <vector>

using namespace std;

/*
 *double funcF(valarray<double>& target, const vector<double>& param);
 *double funcG(valarray<double>& target, const vector<double>& param);
 *double dfuncF(valarray<double>& target, const vector<double>& param);
 *double dfuncG(valarray<double>& target, const vector<double>& param);
 *
 */
double Initial_Condition(valarray<double>& target, const vector<double>& param);

#endif
