#ifndef H_CLSP_H
#define H_CLSP_H

#include <blitz/array.h>
using namespace blitz;

class HeurClsp
{
private:
    int verbose;//verbosity parameter
    int period;//number of period
    int product;//number of product
    int cycle;//maximum number of loop
    double eps;//minimal difference between lower and upper bound
    double param;//smoothing paramter
    //demand function is linear
    Array<double,2>* alpha;//slope
    Array<double,2>* beta;//intercept
    Array<double,2>* prod;//production cost
    Array<double,2>* stor;//holding cost
    Array<double,2>* cons;//consumption of ressouce 
    Array<double,2>* setup;//setup structure
    Array<double,2>* price;//setup structure
    Array<double,2>* production;//setup structure
    Array<double,2>* storage;//setup structure
    Array<int,2>* ind;//index structure from thomas algorithm
    Array<double,1>* constraint;//prodcution constraint
    Array<double,1>* coef;//Khun Thucker coeficient

public:
    //default constructor
    HeurClsp(double* alpha, double* beta, double* prod, double* stor,
        double* consumption, double* setup, double* constraint, int period,
        int product, int verbose, int cycle, double eps, double param);
    double heursolver(); //solver
    void thomas();//dynamic programming solver based on Thomas's paper
    void coefheur();//heuristic who update KKT coef
    double objective();//compute objective
};
#endif
