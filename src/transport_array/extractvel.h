#ifndef EXTRACT_VEL_H_
#define EXTRACT_VEL_H_

#include "advectiveflux.h"
#include "transfunc.h"

// Create a constant test velocity field
int constVelField(vector<vertex>& velocityField, int M, int N);

// Print the velocity field to a txt file
int printVelField(const vector<vertex>& velocityField, int M, int N);

// Extract velocity on quadrature points on a vertical edge
int getQuadVelVert(vertexSet& quadvel, 
                   const vector<vertex>& velocityField, int i, int j, int M, int N);


// Extract velocity on quadrature points on a horizontal edge
int getQuadVelHori(vertexSet& quadVel, 
                   const vector<vertex>& velocityField, int i, int j, int M, int N);

#endif
