//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "QPSolver.hpp"
#include <blitz/array.h>

using namespace Ipopt;
using namespace blitz;

// constructor
QPSolver::QPSolver(Array<double,2> _alpha, Array<double,2> _beta, Array<double,2> _prod,
    Array<double,2> _stor, Array<double,2> _consumption, Array<double,2> _setup,
    Array<double,1> _constraint, int _period, int _product)
{
    period = _period;
    product = _product;
    alpha = new Array<double,2>(_alpha);
    beta = new Array<double,2>(_beta);
    prod = new Array<double,2>(_prod);
    stor = new Array<double,2>(_stor);
    consumption = new Array<double,2>(_consumption);
    setup =  new Array<double, 2>(_setup);
    constraint = new Array<double,1>(_constraint);
    //
    varcoef = new Array<double,1>(period);
    varprod = new Array<double,2>(product,period);
    varstor = new Array<double,2>(product,period);
    varprice = new Array<double,2>(product,period);
}

//destructor
QPSolver::~QPSolver(){}

// returns the size of the problem
bool QPSolver::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style) 
{
    n = 3*period*product; //number of variable
    m = (product+1)*period; //number of constraint
    nnz_jac_g = 5*period*product-1;//number of non zero value in constraint jacobian
    nnz_h_lag = product*period;//number of non zero value in lagrangian hessian
    index_style = TNLP::C_STYLE;// use the C style indexing (0-based)
    return true;
}

// returns the variable bounds
bool QPSolver::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u) 
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    assert(n == 3*period*product);
    assert(m == (product+1)*period);
    for (Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
            x_l[t+j*period] = 0.0; // the variables are positives
            x_u[t+j*period] = sum((*alpha)(j,Range::Range(t,period)))*(*setup)(j,t);  // null in function of setup structure
            x_l[t+j*period+product*period] = 0.0;
            x_u[t+j*period+product*period] = (*alpha)(j,t)/(*beta)(j,t);
            x_l[t+j*period+2*product*period] = 0.0;
            x_u[t+j*period+2*product*period] = 2e19;
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
        g_l[t+period*product] = 0.0;
        g_u[t+period*product] = (*constraint)(t);
    }
    return true;
}

// returns the initial point for the problem
bool QPSolver::get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda)
{
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
    //find the minimal of consumption per time period
    for ( Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
            // initialize to the given starting point
            x[t+j*period] = (*constraint)(t)/
                (max((*consumption)(Range::all(),t))*product); //production minimize the constraint
            //price to stay in feasible region
            x[t+j*period+period*product]=((*alpha)(j,t)-x[t+j*period])/(*beta)(j,t);
            x[t+j*period+2*period*product]=0.;//no storage
        }
    }
    return true;
}

// returns the value of the productective function
bool QPSolver::eval_f(Index n, const Number* x, bool new_x, Number& product_value) 
{
    assert(n == 3*period*product);
    product_value = 0.0;
    for(Index j = 0; j < product; j ++) 
    {
        for (Index t = 0; t < period; t++)
        {
            product_value -= ((*alpha)(j,t)-(*beta)(j,t)*x[period*product+t+j*period])
                *x[period*product+t+j*period]
                -(*prod)(j,t)*x[t+j*period]
                -(*stor)(j,t)*x[2*period*product+t+j*period];
        }
    }
    return true;
}

// return the gradient of the productective function grad_{x} f(x)
bool QPSolver::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) 
{
    assert(n == 3*period*product);
    for(Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
        grad_f[t+j*period] = (*prod)(j,t);
        grad_f[t+j*period+period*product] = 2*(*beta)(j,t)*x[period*product+t+j*period]-(*alpha)(j,t);
        grad_f[t+j*period+2*period*product] = (*stor)(j,t);
        }
    }
    return true;
}

// return the value of the constraints: g(x)
bool QPSolver::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    assert(n == 3*period*product);
    assert(m == (product+1)*period);
    //equatlity constraint
    for(Index j = 0; j < product; j ++)
    {
        for (Index t = 0; t < period; t ++)
        {
            if (t > 0)
            {
                g[t+j*period] = x[t+j*period]
                    + x[2*period*product+t+j*period-1]
                    - x[2*period*product+t+j*period]
                    + (*beta)(j,t)*x[period*product+t+j*period];
            }
            else
            {
                g[t+j*period] =
                    x[t+j*period]
                    - x[2*period*product+t+j*period]
                    + (*beta)(j,t)*x[period*product+t+j*period];
            }
        }
    }
    // inequality constraint
    for (Index t = 0; t < period; t ++)
    {
        g[t+period*product] = 0;
        //printf("\ncosntranit inequality %d\n",t);
        for (Index j = 0; j < product; j ++)
        {
            g[t+period*product] += (*consumption)(j,t)*x[t+j*period];
        }
    }
    return true;
}

// return the structure or values of the jacobian
bool QPSolver::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index *jCol,
                            Number* values) 
{
    //printf("\n");
    if (values == NULL) 
    {
    // return the structure of the jacobian
        //element on 3 diagonals
        for(Index k=0;k<3;k++)
        {
            for(Index t=0; t<period*product;t++) 
            {
                iRow[t+k*period*product]=t;
                jCol[t+k*period*product]=t+k*period*product;
            }
        }
        // element due to the inequality constraint
        for(Index j = 0; j < product; j ++)
        {
            for (Index t = 0; t < period; t ++)
            {
                iRow[t+j*period+3*period*product]=t+period*product;
                jCol[t+j*period+3*period*product]=t+j*period;
            }
        }
        // elements on the sub-diagonale (t_{t-1})
        for(Index t = 0; t < period*product-1; t ++)
        {
            iRow[t+4*period*product]=t+1;
            jCol[t+4*period*product]=t+2*period*product;
        }
    }
    else 
    {
    // return the values of the jacobian of the constraints
        for(Index j = 0; j < product; j ++)
        {
            for (Index t = 0; t < period; t ++)
            {
                values[t+j*period] = 1.;
                values[t+j*period+period*product] = (*beta)(j,t);
                values[t+j*period+2*period*product] = -1;
            }
        }
        // element due to the inequality constraint
        for(Index j = 0; j < product; j ++)
        {
            for (Index t = 0; t < period; t ++)
            {
                values[t+j*period+3*period*product]=(*consumption)(j,t);
            }
        }
        // elements on the sub-diagonale (t_{t-1})
        for(Index t = 0; t < period*product-1; t ++)
        {
            values[t+4*period*product]=1.;
        }
    }
    return true;
}

//return the structure or values of the hessian
bool QPSolver::eval_h(Index n, const Number* x, bool new_x,
                   Number product_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
    if (values == NULL)
    {
        // the hessian for this problem is actually dense
        for (Index t = 0; t <period*product; t ++)
        {
            iRow[t] = t+period*product;
            jCol[t] = t+period*product;
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

void QPSolver::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number product_value,
              const IpoptData* ip_data,
              IpoptCalculatedQuantities* ip_cq) 
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
