#ifndef BRMIXED_H_
#define BRMIXED_H_

class BRMixed {

    public:

        BRMixed();
        ~BRMixed();

        //! Mapping from local degree of freedom to
        //! global index of degree of freedom
        vertex LocalToGlobal(); 


        //! Mapping from global degree of freedom to
        //! local index of degree of freedom
        vertex GlobalToLocal();

    private:



}

#endif
