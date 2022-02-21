rm -rf build
mkdir build
cd build
cmake ..
make
./test -M 17 -N 17
./test -M 33 -N 33
./test -M 65 -N 65
./test -M 129 -N 129
./test -M 257 -N 257
