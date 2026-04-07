cd build
#make; rm *.dat; ./driver -M 50 -N 50 -dt 0.001 -tmax 100 -maxIter 200 -tol 1e-14
make; rm *.dat; ./driver -M 20 -N 20 -dt 0.0001 -tmax 1000 -maxIter 200 -tol 1e-12
