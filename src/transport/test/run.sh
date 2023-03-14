rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Transport Test Built!"
echo " "

valgrind --leak-check=full --show-leak-kinds=all ./test
