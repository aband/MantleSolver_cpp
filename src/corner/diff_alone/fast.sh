cd build
#make; rm *.dat; valgrind ./driver -M 5 -N 5 -dt 1 -tmax 1 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-15
#make; rm *.dat; ./driver -M 20 -N 20 -dt 100 -tmax 2000 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-15
make; rm *.dat; ./driver -M 20 -N 20 -dt 0.1 -tmax 1 -ystart -1.0 -H 2 -xstart -1.0 -L 2 -maxIter 200 -tol 1e-15
