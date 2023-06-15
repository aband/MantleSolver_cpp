#include "util.h"

int main(){

    //derivative d1 = {{1,3.4},{2,5.6},{100,11}};
    derivative d1;
    derivative d2 = {{1,-10},{50,2}};

    cout << "Original maps ..." << endl;

    unordered_map_print(d1);
    unordered_map_print(d2);

    double modi = 0.5;

    unordered_map_arithmetic(d1,d2,std::minus<double>(),modi,std::multiplies<double>());

    cout << "After arithmetic ..." << endl;

    unordered_map_print(d1);
    unordered_map_print(d2);

    return 0;
}
