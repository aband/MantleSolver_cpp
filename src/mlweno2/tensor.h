#ifndef TENSOR_H_
#define TENSOR_H_

/**!
 * Define a tensor object class
 */

#include <cassert>
#include "util.h"

template <class T>
class Tensor{
    public:
        Tensor(const int& r) {setRank(r);};
        ~Tensor() {};

        int setRank(const int& r)
        {rank = r;
         dim.resize(r); 
         val.resize(r); return 1;}

        int setSize(vector<int> inputdim)
        {assert(inputdim.size() == rank);
         dim = inputdim; return 1;}

        int getIndex(const vector<int>& index) 
        {return flattern(index);}

        // Get corresponding values using rank n index
        T &operator()(const vector<int>& index){
            assert(index.size() == rank);
            return val[flattern(index)];
        } 

        T operator()(const vector<int>& index) const{
            assert(index.size() == rank);
            return val[flattern(index)];
        }

    private:
        int rank;
        int size;

        vector<int> dim;
        vector<T> val;

        int flattern(const vector<int>& index){
            int in = 0;
            double multiplier = 1;
            for (int i=0; i<rank; i++){
                assert(index.at(i) < dim.at(i));
                in += index.at(i) * multiplier;
                multiplier *= dim[i];
            }
            return in;
        } 
};

#endif
