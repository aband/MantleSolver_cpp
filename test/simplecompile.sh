#!/bin/bash

gcc testweno.cpp ../src/weno/weno_basis.cpp ../src/integral.cpp -o testweno -lstdc++ -lm
