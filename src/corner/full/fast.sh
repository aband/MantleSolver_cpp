cd build
make; rm *.dat; ./driver -M 40 -N 20 -dt 1 -tmax 1 -ystart -0.5 -H 0.5 -xstart -0.5 -L 1.0 -maxIter 200 -tol 1e-14
