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
    Array<double,1>* discretprices;//vector of allowed prices
    //variable
    Array<double,2>* setup;//setup structure
    Array<double,2>* price;//price
    Array<double,2>* production;//production
    Array<double,2>* storage;//holding
    Array<int,2>* ind;//index structure from thomas algorithm
    Array<double,1>* coef;//Khun Thucker coeficient
    //gap variable as list
    list gap;
    //this pointors allow to change 
    //algorithm behavior
    double (HeurClsp::*updatekkt) ();
    double (HeurClsp::*cost) (blitz::Array<double, 2>, int, int, int);
    double (HeurClsp::*dpprice) (int, int, int);

    double coefheur();//heuristic who update KKT coef
    double coefQP();//QP solver who update KKT coef
    void subproblem();//heurcoef for discret price
    void initVariables();//initiate all variables

    //different specification of price and cost computation for
    //thomas algorithm
    inline double tcost(Array<double,2> tprice, int t, int t0, int j);
    inline double wwcost(Array<double,2> tprice, int t, int t0, int j);
    inline double tprice(int t, int t0, int j);
    inline double wwprice(int t, int t0, int j);
    inline double dprice(int t, int t0, int j);

    //function to interface python and C++ objects
    list ArrayToList(Array<double,2> array);//function to convert blitz array to python list
    void ListToDiscretprices(list prices);//function to import discret prices
public:
    //default constructor
    HeurClsp(list alpha, list beta, list prod, list stor,
        list consumption, list setup, list constraint, int period,
        int product, int verbose, int cycle, float eps, float param);
    //copy constructor
    HeurClsp(const HeurClsp& origin);
    double heursolver();//PCLSP solver
    void thomas();//CLSP solver based on Thomas's paper
    void thomas(list prices);//thomas algorithm with discret prices
    void thomas(double price);//wagner and within algorithm
    double objective();//compute objective
    bool feasible();//return true if the current state are feasible
    void setHeur();//use heurcoef in heursolver

    //methods to get variables;
    list getPrice();
    list getProd();
    list getHold();
    list getSetup();
    list getCoef();
    list getGAP();
    //methods to get and set parameters
    double getSmooth();
    double getStopDiff();
    double getNbCycle();
    int getVerbosity();
    void setSmooth(double param);
    void setStopDiff(double eps);
    void setNbCycle(double cycle);
    void setVerbosity(int verbose);
    //methods to get state
    bool isDiscret();//use discret price
    bool isWW();//use Wagner and within algo

    void plotParam();//plot all parameters
    void plotVariables();//plot all varaibles
};
#endif
