cd build
make; rm *.dat; ./convTest -M 2 -N 80 -dt 1 -tmax 1 -ystart -2 -H 4 -xstart -0.2 -L 0.4 -maxIter 500 -tol 1e-14
#make; rm *.dat; ./convTest -M 2 -N 100 -dt 1 -tmax 1 -ystart -0.25 -H 0.5 -xstart -0.2 -L 0.4
