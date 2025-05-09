cd build
rm *.dat
make

echo " "
echo -n "Start running."
echo " "

./test -M 5 -N 20 -dt 10 -tmax 200

cp CD20.dat restartCD.dat
cp HD20.dat restartHD.dat

echo " "
echo -n "Restart running with 200 time step."
echo " "

./restart -M 5 -N 20 -dt 10 -tmax 200 -start 20

cp CD40.dat restartCD.dat
cp HD40.dat restartHD.dat

echo " "
echo -n "Restart running with 400 time step."
echo " "

./restart -M 5 -N 20 -dt 10 -tmax 200 -start 40

cp CD60.dat restartCD.dat
cp HD60.dat restartHD.dat

echo " "
echo -n "Restart running with 600 time step."
echo " "

./restart -M 5 -N 20 -dt 10 -tmax 200 -start 60
