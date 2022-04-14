rm -rf build
mkdir build
cd build
cmake ..
make
#./test -meshtype 1 -M 2 -N 2 -printmesh 0
#./test -meshtype 0 -M 16 -N 16
#./test -meshtype 0 -M 32 -N 32
./test -meshtype 0 -M 10 -N 10
#./test -meshtype 0 -M 128 -N 128
#./test -meshtype 0 -M 256 -N 256
#./test -meshtype 0 -M 512 -N 512
