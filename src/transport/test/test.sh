rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Transport Test Built!"
echo " "

#valgrind --leak-check=full --show-leak-kinds=all ./test

# A half decent run
./test -M 10 -N 10 -Tmax 0.001 -dt 0.001 -implicit 1 -snes_monitor
