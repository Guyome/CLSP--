#ifndef H_CLSP_H
#define H_CLSP_H

#include <boost/python.hpp>
#include <blitz/array.h>
using namespace blitz;
using namespace boost::python;

class HeurClsp
{
private:
    //parameters
    int verbose;//verbosity parameter
    int period;//number of period
    int product;//number of product
    int cycle;//maximum number of loop
    double eps;//minimal difference between lower and upper bound
    double param;//smoothing paramter
    //data
    //demand function is linear
    Array<double,2>* alpha;//slope
    Array<double,2>* beta;//intercept
    Array<double,2>* prod;//production cost
    Array<double,2>* stor;//holding cost
    Array<double,2>* cons;//consumption of ressouce 
    Array<double,2>* setupcost;//setup cost
    Array<double,1>* constraint;//prodcution constraint
    //variable
    Array<double,2>* setup;//setup structure
    Array<double,2>* price;//price
    Array<double,2>* production;//production
    Array<double,2>* storage;//holding
    Array<int,2>* ind;//index structure from thomas algorithm
    Array<double,1>* coef;//Khun Thucker coeficient

    void (HeurClsp::*updatekkt) ();
    double (HeurClsp::*cost) (blitz::Array<double, 2>, int, int, int);
    double (HeurClsp::*dpprice) (int, int, int);

    void coefheur();//heuristic who update KKT coef
    void coefQP();//QP solver who update KKT coef
    void subproblem();//heurcoef for discret price
    void initVariables();//initiate all variables
    double tcost(Array<double,2> tprice, int t, int t0, int j);
    double wwcost(Array<double,2> tprice, int t, int t0, int j);
    double tprice(int t, int t0, int j);
    double wwprice(int t, int t0, int j);

    list ArrayToList(Array<double,2> array);//function to convert blitz array to python list

public:
    //default constructor
    HeurClsp(list alpha, list beta, list prod, list stor,
        list consumption, list setup, list constraint, int period,
        int product, int verbose, int cycle, float eps, float param);
    //copy constructor
    HeurClsp(const HeurClsp& origin);
    double heursolver();//PCLSP solver
    void thomas();//CLSP solver based on Thomas's paper
    double objective();//compute objective
    bool feasible();//return true if the current state are feasible
    void setHeur();//use heurcoef in heursolver
    double ww(double price);//wagner and within algorithm

    //methods to get variables;
    list getPrice();
    list getProd();
    list getHold();
    list getSetup();
    list getCoef();

    void plotParam();//plot all parameters
    void plotVariables();//plot all varaibles
};
#endif
