#ifndef __QPSolver_HPP__
#define __QPSolver_HPP__

#include <coin/IpTNLP.hpp>
#include <blitz/array.h>

using namespace Ipopt;
using namespace blitz;

/**
 * @brief Specification of Ipopt abstract class (see https://projects.coin-or.org/Ipopt)
 * for solve QP program due to KKT coefficient update in HeurClsp()
 *
 * @author Guillaume Lanquepin <guillaume@himolde.no>
 * @version 0.1
 */
class QPSolver : public TNLP
{
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
    * @param product the number of goods*/    
    QPSolver(Array<double,2>* alpha, Array<double,2>* beta, Array<double,2>* prod, Array<double,2>* stor,
                    Array<double,2>* consumption, Array<double,2>* setup, Array<double,1>* constraint,
                    int period, int product);
    virtual ~QPSolver();/**< Default destructor */

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                    Index& nnz_h_lag, IndexStyleEnum& index_style);
    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                    Index m, Number* g_l, Number* g_u);
    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda);
    /** Method to return the objective value */
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
    /** Method to return the constraint residuals */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
    /** Method to return:
    *   - The structure of the jacobian (if "values" is NULL)
    *   - The values of the jacobian (if "values" is not NULL)
    */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                                Index m, Index nele_jac, Index* iRow, Index *jCol,
                                Number* values);
    /** Method to return:
    *   - The structure of the hessian of the lagrangian (if "values" is NULL)
    *   - The values of the hessian of the lagrangian (if "values" is not NULL)
    */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                                Number obj_factor, Index m, const Number* lambda,
                                bool new_lambda, Index nele_hess, Index* iRow,
                                Index* jCol, Number* values);
    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(SolverReturn status,
                                Index n, const Number* x, const Number* z_L, const Number* z_U,
                                Index m, const Number* g, const Number* lambda,
                                Number obj_value,const IpoptData* ip_data,
                                IpoptCalculatedQuantities* ip_cq);
    //@}
    
    /** @name Attributes Methods
    * Following functions are initiate by finalize_solution(), thus they return NULL
    * before optimization ( see OptimizeTNLP() function in IPOPT)*/
    //@{
    Array<double,1> getCoef();
    Array<double,2> getPrice();
    Array<double,2> getProd();
    Array<double,2> getStor();
    //@}

    private:
    /**@name Methods to block default compiler methods.
    * The compiler automatically generates the following three methods.
    *  Since the default compiler implementation is generally not what
    *  you want (for all but the most simple classes), we usually
    *  put the declarations of these methods in the private section
    *  and never implement them. This prevents the compiler from
    *  implementing an incorrect "default" behavior without us
    *  knowing. (See Scott Meyers book, "Effective C++")
    *
    */
    //@{
    /** QPSolver(); */
    QPSolver(const QPSolver&);
    QPSolver& operator=(const QPSolver&);
    //@}
    
    /**@name QP parameters and variables*/
    //@{
    int period;/**< Number of time period */
    int product;/**< Number of product */
    Array<double,2>* alpha;/**< Slope of demand function */
    Array<double,2>* beta;/**< Intercept of demand function */
    Array<double,2>* prod;/**< Production cost */
    Array<double,2>* stor;/**< Holding cost */
    Array<double,2>* consumption;/**< Consumption of resource per product*/
    Array<double,2>* setup;/**< Setup structure */
    Array<double,1>* constraint;/**< Production constraint */
    Array<double,1>* varcoef;/**< Khun Thucker coefficient */
    Array<double,2>* varprice;/**< Optimal price */
    Array<double,2>* varprod;/**< Optimal production */
    Array<double,2>* varstor;/**< Optimal storage */
    //@}
};


#endif
