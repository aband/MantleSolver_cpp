#ifndef EXTRACT_VEL_H_
#define EXTRACT_VEL_H_

// Create a constant test velocity field
int constVelField(vector<vertex>& velocityField, int M, int N);

// Extract velocity on quadrature points on a vertical edge
int getQuadVelVert(vertexSet& quadvel, 
                   const vector<vertex>& velocityField, int i, int j, int M, int N);


// Extract velocity on quadrature points on a horizontal edge
int getQuadVelHori(vertexSet& quadVel, 
                   const vector<vertex>& velocityField, int i, int j, int M, int N);

#endif
