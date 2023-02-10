rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Phase test case Built!"
echo " "

#valgrind --leak-check=full ./test
