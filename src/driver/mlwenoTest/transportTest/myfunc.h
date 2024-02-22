#ifndef MYFUNC_H_
#define MYFUNC_H_

enum Location {leftBndry, rightBndry, topBndry, bottomBndry, interior};

const std::array<double,2> advFunc(double u);

const std::array<double,2> dAdvFunc(double u);

const double diffFunc(double u);

const double dDiffFunc(double u);

double Distribution(const vertex& point,
                    const vector<double>& param);

// Assign location to mlwenouse object 
Location assignLocation(const indice& globalCell);

#endif
