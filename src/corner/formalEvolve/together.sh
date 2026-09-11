cd build
#make; rm *.dat; ./preheat -M 20 -N 20 -dt 25 -tmax 1000 -frame 10 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 200 -tol 1e-14
#cp /home/renpo/Research/MantleSolver_cpp/src/corner/formalpreheat/archiveSol/20/* . 
cp /home/renpo/Research/MantleSolver_cpp/src/corner/formalpreheat/build/* . 
make; ./together -M 32 -N 32 -dt 5 -tmax 50 -frame 5 -ystart -0.5 -H 0.5 -xstart 0.0 -L 0.5 -maxIter 100 -tol 1e-14
