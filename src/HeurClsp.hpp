#ifndef H_CLSP_H
#define H_CLSP_H

#include <blitz/array.h>
using namespace blitz;

class HeurClsp
{
private:
    int verbose//verbosity parameter
    int period;//number of period
    int product;//number of product
    //demand function is linear
    Array<double,2> alpha;//slope
    Array<double,2> beta;//intercept
    Array<double,2> prod;//production cost
    Array<double,2> stor;//holding cost
    Array<double,2> cons;//consumption of ressouce 
    Array<double,2> setup;//setup structure
    Array<double,1> constraint;//prodcution constraint
    Array<double,1> coef;//Khun Thucker coeficient

protected:
    //dynamic programming solver based on Thomas's paper
    Array<int,2> thomas(Array<double,2> results);
    
public:
    //default constructor
    HeurClsp(Array<double,2> alpha, Array<double,2> beta, Array<double,2> prod, Array<double,2> stor,
        Array<double,2> consumption, Array<double,2> setup, Array<double,1> constraint, int period,
        int product);
}
#endif
