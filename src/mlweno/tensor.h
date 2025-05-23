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
        Tensor () {};
        Tensor(const int& r) {setRank(r);};
        ~Tensor() {};

        int setRank(const int& r)
        {rank = r;
         dim.resize(r); return 1;}

        int setSize(vector<int> inputdim)
        {//assert(inputdim.size() == rank);
         rank = inputdim.size();
         dim = inputdim; 
         size = 1; for (const auto& it : dim) {size *= it;}
         val.resize(size);
         return 1;}

        int getIndex(const vector<int>& index) 
        {return flattern(index);}

        int getSize() const
        {return size;}

        int getSize(const int& i) const
        {return dim.at(i);}

        // Get corresponding values using rank n index
        T &operator()(const vector<int>& index){
            assert(index.size() == rank);
            return val[flattern(index)];
        } 

        T operator()(const vector<int>& index) const{
            assert(index.size() == rank);
            return val[flattern(index)];
        }

        T &operator()(const int& index){
            assert(index < size);
            return val[index];
        } 

        T operator()(const int& index) const{
            assert(index < size);
            return val[index];
        }

    private:
        int rank;
        int size;

        vector<int> dim;

        int flattern(const vector<int>& index) const{
            int in = 0;
            int multiplier = 1;
            for (int i=0; i<rank; i++){
                assert(index.at(i) < dim.at(i));
                in += index.at(i) * multiplier;
                multiplier *= dim[i];
            }
            return in;
        }

        vector<T> val;
};

// Some arithmetic functions that will be used
// Restricted to double data type for the fact that T can be 
// non numeric type like pointer
int Tensor_add(const Tensor<double>& t1, 
               const Tensor<double>& t2, 
               Tensor<double>& t3);

int Tensor_multi(const Tensor<double>& t1, 
                 const Tensor<double>& t2, 
                 Tensor<double>& t3);

int Tensor_zero(Tensor<double>& t);

int Tensor_multi_add(const Tensor<double>& t1, 
                     const Tensor<double>& t2, 
                     double scale,
                     Tensor<double>& t3);

int Tensor_scale(double scale, Tensor<double>& t);

#endif
