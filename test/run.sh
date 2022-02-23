rm -rf build
mkdir build
cd build
cmake ..
make
#./test -meshtype 1 -M 2 -N 2 -printmesh 0
./test -meshtype 1 -M 17 -N 17
./test -meshtype 1 -M 33 -N 33
./test -meshtype 1 -M 65 -N 65
./test -meshtype 1 -M 129 -N 129
./test -meshtype 1 -M 257 -N 257
./test -meshtype 1 -M 513 -N 513
