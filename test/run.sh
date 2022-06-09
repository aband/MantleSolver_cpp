rm -rf build
mkdir build
cd build
cmake ..
make

echo " "
echo -n "Multilevel Weno"
echo " "
./multilevel_weno -M 5 -N 5
#./multilevel_weno -M 40 -N 40
#./multilevel_weno -M 80 -N 80
#./multilevel_weno -M 160 -N 160

#echo " "
#echo -n "Original Weno"
#echo " "
#./original_weno -meshtype 0 -M 20 -N 20
#./original_weno -meshtype 0 -M 40 -N 40

#./original_weno
#./original_weno -meshtype 0 -M 60 -N 60
#./original_weno -meshtype 0 -M 150 -N 150
#./original_weno -meshtype 0 -M 129 -N 129
#./original_weno -meshtype 1 -M 257 -N 257
#./original_weno -meshtype 1 -M 513 -N 513
