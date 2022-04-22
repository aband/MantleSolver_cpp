rm -rf build
mkdir build
cd build
cmake ..
make
#./test -meshtype 1 -M 2 -N 2 -printmesh 0
#./test -meshtype 0 -M 7 -N 7
./test -meshtype 0 -M 60 -N 60
#./test -meshtype 0 -M 150 -N 150
#./test -meshtype 0 -M 129 -N 129
#./test -meshtype 1 -M 257 -N 257
#./test -meshtype 1 -M 513 -N 513
