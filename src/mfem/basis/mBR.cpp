#include <../../include/FEbasis.h>

/*
 *Define Modified Bernardi-Raugel basis function
 */

double r(Element& elem, point& x)
{
    Vector2D x0, x1, x2, x3;
           
    x0 = (*elem).corner[0];
    x1 = (*elem).corner[1];
    x2 = (*elem).corner[2];
    x3 = (*elem).corner[3];

    return lambda(elem,x,2)*lambda(elem,x,3)/(lambda(elem,x0,2)*lambda(elem,x0,3)) - 
           lambda(elem,x,3)*lambda(elem,x,0)/(lambda(elem,x1,3)*lambda(elem,x1,0)) +
           lambda(elem,x,0)*lambda(elem,x,1)/(lambda(elem,x2,0)*lambda(elem,x2,1)) -
           lambda(elem,x,1)*lambda(elem,x,2)/(lambda(elem,x3,1)*lambda(elem,x3,2)) ;
}


