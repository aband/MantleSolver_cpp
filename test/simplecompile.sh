#!/bin/bash

gcc testweno.cpp ../weno/weno_basis.cpp ../integral/integral.cpp -o testweno -lstdc++ -lm
