//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "QPSolver.hpp"
#include <blitz/array.h>

using namespace Ipopt;
using namespace blitz;

//constructor
QPSolver::QPSolver(Array<double,2>* _alpha, Array<double,2>* _beta, Array<double,2>* _prod,
    Array<double,2>* _stor, Array<double,2>* _consumption, Array<double,2> *_setup,
    Array<double,1>* _constraint, int _period, int _product)
{
    //problem parameters
    period = _period;
    product = _product;
    alpha = _alpha;
    beta = _beta;
    prod = _prod;
    stor = _stor;
    consumption = _consumption;
    setup = _setup;
    constraint = _constraint;
    //only needed to get back results
    varcoef = new Array<double,1>(period);
    varprod = new Array<double,2>(product,period);
    varstor = new Array<double,2>(product,period);
    varprice = new Array<double,2>(product,period);
}

//destructor
QPSolver::~QPSolver(){}

//returns the size of the problem
bool QPSolver::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
    Index& nnz_h_lag, IndexStyleEnum& index_style) 
{
    n = 3*period*product;//number of variable
    m = (product+1)*period;//number of constraint
    nnz_jac_g = 5*period*product-1;//number of non zero value in constraint jacobian
    nnz_h_lag = product*period;//number of non zero value in lagrangian hessian
    index_style = TNLP::C_STYLE;//use the C style indexing (0-based)
    return true;
}

//returns the variable bounds
bool QPSolver::get_bounds_info(Index n, Number* x_l, Number* x_u,
    Index m, Number* g_l, Number* g_u) 
{
    for (Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
            //the variables are positive
            x_l[t+j*period] = 0;
            x_l[t+j*period+product*period] = 0;
            x_l[t+j*period+2*product*period] = 0;
            x_u[t+j*period] = 1e19*(*setup)(j,t);//null in function of setup structure
            x_u[t+j*period+product*period] = (*alpha)(j,t)/(*beta)(j,t);//limits price to have a least no demand
            x_u[t+j*period+2*product*period] = 1e19;
        }
    }

    // product*period equality constraint
    for (Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
            g_l[t+j*period] = g_u[t+j*period] = (*alpha)(j,t);
        }
    }
    //inequality constraints
    for (Index t=0; t<period; t++)
    {
        g_l[t+period*product] = 0;
        g_u[t+period*product] = (*constraint)(t);
    }
    return true;
}

//returns the initial point for the problem
bool QPSolver::get_starting_point(Index n, bool init_x, Number* x,
    bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda,
    Number* lambda)
{
    assert(init_x == true);//primal
    assert(init_z == false);//no dual starting point
    assert(init_lambda == false);//no kkt coef initials values
    //find the minimal of consumption per time period
    for ( Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
            //initialize to the given starting point
            x[t+j*period] = ((*constraint)(t)/
                (max((*consumption)(Range::all(),t))*product))*(*setup)(j,t);//production minimize the constraint
            //price to stay in feasible region
            x[t+j*period+period*product]=((*alpha)(j,t)-x[t+j*period])/(*beta)(j,t);
            x[t+j*period+2*period*product]=0;//no storage
        }
    }
    return true;
}

//returns the value of the productective function
bool QPSolver::eval_f(Index n, const Number* x, bool new_x, Number& product_value) 
{
    product_value = 0;
    for(Index j = 0; j < product; j ++) 
    {
        for (Index t = 0; t < period; t++)
        {
            product_value -= ((*alpha)(j,t)-(*beta)(j,t)*x[t+j*period+period*product])*x[t+j*period+period*product]
                -(*prod)(j,t)*x[t+j*period]-(*stor)(j,t)*x[t+j*period+2*period*product];
        }
    }
    return true;
}

//return the gradient of the productective function grad_{x} f(x)
bool QPSolver::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) 
{
    for(Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
        grad_f[t+j*period] = (*prod)(j,t);
        grad_f[t+j*period+period*product] = 2*(*beta)(j,t)*x[t+j*period+period*product]-(*alpha)(j,t);
        grad_f[t+j*period+2*period*product] = (*stor)(j,t);
        }
    }
    return true;
}

//return the value of the constraints: g(x)
bool QPSolver::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    //equatlity constraint
    for(Index j = 0; j < product; j ++)
    {
        g[j*period] = x[j*period] - x[j*period+2*period*product]
                + (*beta)(j,0)*x[j*period+period*product];
        for (Index t = 1; t < period; t ++)
        {
            g[t+j*period] = x[t+j*period] - x[t+j*period+2*period*product]
                + (*beta)(j,t)*x[t+j*period+period*product]
                + x[t-1+j*period+2*period*product];
        }
    }
    // inequality constraint
    for (Index t = 0; t < period; t ++)
    {
        g[t+period*product] = 0;
        for (Index j = 0; j < product; j ++)
        {
            g[t+period*product] += (*consumption)(j,t)*x[t+j*period];
        }
    }
    return true;
}

//return the structure or values of the jacobian
bool QPSolver::eval_jac_g(Index n, const Number* x, bool new_x,
    Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values) 
{
    if (values == NULL) 
    {
        //return the structure of the jacobian
        //elements on 3 diagonals
        for(Index j = 0; j < product; j ++)
        {
            for(Index t = 0; t < period; t++) 
            {
                //derivate by production
                iRow[t+j*period] = t+j*period;
                jCol[t+j*period] = t+j*period;
                iRow[t+j*period+3*period*product] = t+period*product;//ineguality constraint
                jCol[t+j*period+3*period*product] = t+j*period;//ineguality constraint
                //derivate by price
                iRow[t+j*period+period*product] = t+j*period;
                jCol[t+j*period+period*product] = t+j*period+period*product;
                //derivate by storage
                iRow[t+j*period+2*period*product] = t+j*period;
                jCol[t+j*period+2*period*product] = t+j*period+2*period*product;
                if (t+j*period < period*product-1)
                {
                    iRow[t+j*period+4*period*product] = t+j*period+1;
                    jCol[t+j*period+4*period*product] = t+j*period+2*period*product;
                }
            }
        }
    }
    else 
    {
        //return the values of the jacobian of the constraints
        for(Index j = 0; j < product; j ++)
        {
            for (Index t = 0; t < period; t ++)
            {
                //derivate by production
                values[t+j*period] = 1;
                values[t+j*period+3*period*product]= (*consumption)(j,t);//ineguality constraint
                //derivate by price
                values[t+j*period+period*product] = (*beta)(j,t);
                //derivate by storage
                values[t+j*period+2*period*product] = -1;
                if (t+j*period < period*product-1)
                {
                    values[t+j*period+4*period*product] = 1;
                }
            }
        }
    }
    return true;
}

//return the structure or values of the hessian
bool QPSolver::eval_h(Index n, const Number* x, bool new_x,
    Number product_factor, Index m, const Number* lambda,
    bool new_lambda, Index nele_hess, Index* iRow,Index* jCol,
    Number* values)
{
    if (values == NULL)
    {
        //the hessian for this problem is actually dense
        for (Index j = 0; j < product; j ++)
        {
            for (Index t = 0; t < period; t ++)
            {
                iRow[t+j*period]=t+j*period+period*product;
                jCol[t+j*period]=t+j*period+period*product;
            }
        }
    }
    else
    {
        for (Index j = 0; j < product; j ++)
        {
            for (Index t = 0; t < period; t ++)
            {
                values[t+j*period] = 2*(*beta)(j,t)*product_factor;
            }
        }
    }
    return true;
}

//return optimal variables
void QPSolver::finalize_solution(SolverReturn status,
    Index n, const Number* x, const Number* z_L, const Number* z_U,
    Index m, const Number* g, const Number* lambda, Number product_value,
    const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq) 
{
    for (int t = 0; t < period; t ++)
    {
        (*varcoef)(t) = lambda[t+period*product];
    }
    for (int j = 0; j < product; j ++)
    {
        for (int t = 0; t < period; t += 1)
        {
            (*varprod)(j,t) = x[t+j*period];
            (*varprice)(j,t) = x[t+j*period+period*product];
            (*varstor)(j,t) = x[t+j*period+2*period*product];
        }
    }
}

//functions to get back variables
Array<double,1> QPSolver::getCoef(){
    return (*varcoef);
}

Array<double,2> QPSolver::getPrice(){
    return (*varprice);
}

Array<double,2> QPSolver::getProd(){
    return (*varprod);
}

Array<double,2> QPSolver::getStor(){
    return (*varstor);
}
