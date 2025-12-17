cd build
make; rm *.dat; ./midRidge -M 3 -N 3 -dt 1 -tmax 1 -ystart -0.5 -H 0.5 -xstart -0.5 -L 1.0 -maxIter 1000 -tol 1e-17
