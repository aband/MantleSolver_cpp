rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Convergence Test Built!"
echo " "

./test
#valgrind --leak-check=full ./test
#mpiexec -n 2 ./test -M 8 -N 8
#./test -M 11 -N 11
#./test -M 21 -N 21
#./test -M 41 -N 41
