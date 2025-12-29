cd build
make; rm *.dat; ./driver -M 5 -N 5 -dt 1 -tmax 1 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-15 -transport 1
