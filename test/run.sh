rm -rf build
mkdir build
cd build
cmake ..
make
#./test -meshtype 1 -M 2 -N 2 -printmesh 0
#./test -meshtype 0 -M 9 -N 9
./test -meshtype 0 -M 16 -N 16
#./test -meshtype 0 -M 65 -N 65
#./test -meshtype 0 -M 129 -N 129
#./test -meshtype 1 -M 257 -N 257
#./test -meshtype 1 -M 513 -N 513
