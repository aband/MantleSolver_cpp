rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

echo " "
echo -n "Convergence Test Built!"
echo " "

#./test -M 20 -N 20
#valgrind --leak-check=full -s ./test
#mpiexec -n 2 ./test -M 8 -N 8
#./test -M 11 -N 11
#./test -M 21 -N 21
#./test -M 41 -N 41
