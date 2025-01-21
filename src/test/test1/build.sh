rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Test Built Successfully!"
echo " "

#valgrind --leak-check=full --show-leak-kinds=all ./test

# A half decent run
#./test -M 50 -N 50 -Tmax 0.5 -dt 0.005
