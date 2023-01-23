rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Convergence Test Built!"
echo " "

valgrind ./test
#./test -M 11 -N 11
#./test -M 21 -N 21
#./test -M 41 -N 41
