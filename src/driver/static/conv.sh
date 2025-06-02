cd build
make
./convTest -M 2 -N 20 -dt 1 -tmax 1 -ystart -2 -H 4 -xstart -0.2 -L 0.4 -maxIter 500 -tol 1e-14
./convTest -M 2 -N 40 -dt 1 -tmax 1 -ystart -2 -H 4 -xstart -0.2 -L 0.4 -maxIter 500 -tol 1e-14
./convTest -M 2 -N 80 -dt 1 -tmax 1 -ystart -2 -H 4 -xstart -0.2 -L 0.4 -maxIter 500 -tol 1e-14
./convTest -M 2 -N 160 -dt 1 -tmax 1 -ystart -2 -H 4 -xstart -0.2 -L 0.4 -maxIter 500 -tol 1e-14



