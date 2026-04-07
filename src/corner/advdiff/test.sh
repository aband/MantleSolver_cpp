cd build
make; rm *.dat; ./driver -M 20 -N 20 -dt 100 -tmax 800 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-14
