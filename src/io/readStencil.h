#ifndef READSTENCIL_H_
#define READSTENCIL_H_

#include <iostream>

#include "reader.h"

// Read stencils from input file
int readWenoStencil(std::vector<std::array<int,2>>,
                    std::vector<std::vector<std::array<int,2>>>);

#endif
