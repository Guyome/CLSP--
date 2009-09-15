#ifndef H_CLSP_H
#define H_CLSP_H

#include <boost/python.hpp>
#include <blitz/array.h>
using namespace blitz;
using namespace boost::python;

/**
 * @brief This library aims to solve lot-sizing related problems,
 * with the following shape for J goods and T as time horizons:
 * \f[ \sum_{j=1}^J\sum_{t=1}^T d_{j,t}p_{j,t}-x_{j,t}v_{j,t}-
 * h_{j,t}i_{j,t}-s_{j,t}\delta_{j,t} \f]
 * s.c. \f[ d_{j,t}=x_{j,t}+i_{j,t-1}-i{j,t}\f]
 *   where:
 *       - \f$ d_{j,t}\f$ is the demand
 *       - \f$ p_{j,t}\f$ is the price
 *       - \f$ x_{j,t}\f$ is the production
 *       - \f$ h_{j,t}\f$ is the holding
 *       - \f$ \delta_{j,t}\f$ is the current
 * setup (integer variable)
 *       - \f$ s_{j,t}\f$,\f$ v_{j,t}\f$,\f$ i_{j,t}\f$ are
 * respectively setup, production and inventory cost and are given
 * It currently manages linear profit maximizing lot-sizing problems
 * (\f$ d_{j,t}=\alpha_{j,t}p_{j,t}-\beta_{j,t}\f$)
 * and capacitated problems (\f$ \sum_{j=1}^J x_{j,t} \leq R_t\f$)
 *
 * @author Guillaume Lanquepin <guillaume@himolde.no>
 * @version 0.1
 */
class HeurClsp
{
private:
    /**@name Problem parameters */
    //@{
    int verbose;/**< verbosity parameter */
    int period;/**< Number of period */
    int product;/**< Number of product */
    int cycle;/**< Maximum number of loop */
    double eps;/**< Minimal difference between lower and upper bound */
    double param;/**< Smoothing parameter */
    //@}


    /**@name Data */
    //@{
    Array<double,2>* alpha;/**< Slope of demand function */
    Array<double,2>* beta;/**< Intercept of demand function */
    Array<double,2>* prod;/**< Production cost */
    Array<double,2>* stor;/**< Holding cost */
    Array<double,2>* cons;/**< Consumption of resource per product*/
    Array<double,2>* setupcost;/**< Setup cost */
    Array<double,1>* constraint;/**< Production constraint per time period*/
    Array<double,1>* discretprices;/**< Vector of allowed prices */
    //@}

    /**@name Variable */
    //@{
    Array<double,2>* setup;/**< Setup structure */
    Array<double,2>* price;/**< Price */
    Array<double,2>* production;/**< Production */
    Array<double,2>* storage;/**< Holding */
    Array<int,2>* ind;/**< Index structure from Thomas's algorithm */
    Array<double,1>* coef;/**< Khun Thucker coefficients */

    list gap;/**< GAP variable as list */
    //@}

    /**@name Thomas() specification
    * This pointers allow to change algorithm behavior by specifying Thomas() 
    * computation for price and holding */
    //@{
    double (HeurClsp::*updatekkt) ();/**< Can be based on heuristic or quadratic problem solver*/
    double (HeurClsp::*cost) (blitz::Array<double, 2>, int, int, int);/**< Can be tcost() or wwcost()*/
    double (HeurClsp::*dpprice) (int, int, int);/**< Can be tprice(),dprice() or wwprice()*/
    double coefheur();/**< Heuristic who update KKT coef */
    double coefQP();/**< QP solver who update KKT coef */
    void subproblem();/**< Heurcoef for discrete price */
    void initVariables();/**< Initiate all variables */
    //@}

    /**@name Thomas price and cost specification
    * Different specification of price and cost computation for
    * Thomas algorithm */
    //@{
    inline double tcost(Array<double,2> tprice, int t, int t0, int j);/**< Return cost for Thomas algo*/
    inline double wwcost(Array<double,2> tprice, int t, int t0, int j);/**< Return cost for Wagner and Within algo.*/
    inline double tprice(int t, int t0, int j);/**< Return price for Thomas algo*/
    inline double wwprice(int t, int t0, int j);/**< Return price for Wagner and Within algo.*/
    inline double dprice(int t, int t0, int j);/**< Return the best price*/
    //@}

    /**@name Python interface
    * Functions to interface python and C++ objects */
    //@{
    list ArrayToList(Array<double,2> array);/**< Function to convert blitz array to python list */
    void ListToDiscretprices(list prices);/**< Function to import discret prices */
    //@}
public:
    /** Default constructor 
    * @param alpha the slope of demand function
    * @param beta the intercept of demand function
    * @param prod the production cost
    * @param stor the holding cost
    * @param consumption the resource consumption per product
    * @param setup the setup cost
    * @param constraint the production constraint per time period
    * @param period the number of time period
    * @param product the number of goods
    * @param verbose the verbosity level (less than 3)
    * @param cycle the maximal number of iteration
    * @param eps the minimal difference between bounds
    * @param param the smoothing coefficient*/
    HeurClsp(list alpha, list beta, list prod, list stor,
        list consumption, list setup, list constraint, int period,
        int product, int verbose, int cycle, float eps, float param);
    HeurClsp(const HeurClsp& origin);/**< Copy constructor 
                                        @param origin ther original HeurClsp object*/
    double heursolver();/**< PCLSP Solver */
    void thomas();/**< CLSP solver based on Thomas's paper */
    void thomas(list prices);/**< CLSP solver algorithm with discrete prices 
                                    @param prices all the allowed prices*/
    void thomas(double price);/**< Wagner and within algorithm
                                    (fixed intercept and null slope as demand function)
                                    @param price the price value in WW algo.*/
    double objective();/**< Return objective */
    bool feasible();/**< Return true if the current state are feasible */
    /** Use heuristic to update KKT coefficient in spite of QP methods */
    void setHeur();

    /** @name Informative methods
    * Methods to obtain variables and various quality information; */
    //@{
    list getPrice();/**< Returns the current price values as list*/
    list getProd();/**< Returns the current production values as list*/
    list getHold();/**< Returns the current holding values as list*/
    list getSetup();/**< Returns the current setup values as list*/
    list getCoef();/**< Returns the current KKT coefficient values as list */
    list getGAP();/**< Returns the current GAP values as list */
    //@}

    /* @name Parametrize methods */
    //@{
	double getProduct();
    double getPeriod();
    double getSmooth();/**< Returns KKT coefficients smoothing parameter*/
    double getStopDiff();/**< Returns bound difference stopping condition */
    double getNbCycle();/**< Returns number of iteration stopping condition*/
    int getVerbosity();/**< Returns the current verbosity*/
    void setSmooth(double param);/**< Set KKT coefficients smoothing parameter
                                        @param param the smoothing coefficient*/
    void setStopDiff(double eps);/**< Set bound difference stopping condition
                                        @param eps the minimal difference between bounds*/
    void setNbCycle(double cycle);/**< Set number of iteration stopping condition
                                        @param cycle the maximal number of iteration*/
    void setVerbosity(int verbose);/**< Set verbosity
                                        @param verbose the verbosity level (less than 3)*/
    //@}

    /** @name State methods
    * methods to get the current state */
    //@{
    bool isDiscret();/**< Use discrete price */
    bool isWW();/**< Use Wagner and within algorithm */
    //@}

    /** @name Plotting methods */
    //@{
    void plotParam();/**< Plot all parameters */
    void plotVariables();/**< Plot all variables */
    //@}
};
#endif
