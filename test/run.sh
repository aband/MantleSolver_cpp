rm -rf build
mkdir build
cd build
cmake ..
make

./multilevel_weno
#./original_weno
#./original_weno -meshtype 0 -M 60 -N 60
#./original_weno -meshtype 0 -M 150 -N 150
#./original_weno -meshtype 0 -M 129 -N 129
#./original_weno -meshtype 1 -M 257 -N 257
#./original_weno -meshtype 1 -M 513 -N 513
