cd build
#make; rm *.dat; ./preheat -M 20 -N 20 -dt 25 -tmax 1000 -frame 10 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-14
./preheat -M 40 -N 40 -dt 10 -tmax 1500 -frame 10 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-14
