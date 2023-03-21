rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Transport Test Built!"
echo " "

#./test

valgrind --leak-check=full --show-leak-kinds=all ./test

# A half decent run
#./test -M 40 -N 40 -Tmax 0.5 -dt 0.01
