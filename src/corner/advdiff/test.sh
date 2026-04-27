cd build
make; rm *.dat; ./driver -M 50 -N 50 -dt 25 -tmax 4000 -frame 50 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-14
