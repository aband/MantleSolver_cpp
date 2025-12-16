cd build
make; rm *.dat; ./midRidge -M 5 -N 5 -dt 1 -tmax 1 -ystart -0.5 -H 0.5 -xstart -0.5 -L 1.0 -maxIter 1000 -tol 1e-17
