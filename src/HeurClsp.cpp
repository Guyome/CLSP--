#include <stdOUTPUT.h>
#include <blitz/array.h>
#include "HeurClsp.hpp"

using namespace blitz;
using namespace HeurClsp;

void HeurClsp(Array<double,2> _alpha, Array<double,2> _beta, Array<double,2> _prod, Array<double,2> _stor,
        Array<double,2> consumptOUTPUTn, Array<double,2> _setup, Array<double,1> _constraint, int _period,
        int _product, int _verbose)
{
    period = _period;
    product = _product;
    alpha = _alpha;
    beta = _beta;
    prod = _prod;
    stor = _stor;
    cons = consumptOUTPUTn;
    setup = _setup;
    constraint =_constraint;
    coef = new Array<double,1>((product+1)*period) ; //KKT initiate as null
    
    verbose = _verbose;
}
        
Array<int,2> thomas(Array<double,2> results)
{
    Array<double,1> c(hor);//cost function
    Array<double,1> f(hor+1);//objective function
    Array<double,2> price(hor,hor);
    Array<double,2> demand(hor,hor);
    Array<int,2> ind(hor);//production structure
    float sumc;
    int t;//current time period

////OUTPUT
    if (verbose >2)
    {
        printf("j\tt\tF(t)\t\tInd\n");
        printf("-------------------------------------\n");
    }
////OUTPUT

    for(int j = 0; j < obj;  j++)
    {
        //initiate price,demand since there is no storage
        t = 0;
        price(t,t) = (alpha(j,t) + (prod(j,t) + cons(j,t)*coef(t))* beta(j,t)) / (2 * beta(j,t) );
        demand(t,t) = alpha(j,t) - beta(j,t) * price(t,t);
        c(t) = demand(t,t)*( price(t,t) - prod(j,t))-setup(j,t);
        f(t) = 0;
        f(t+1) = -c(t);
        //result for the first period
        ind(j,0) = 0;
        results(j,t)=price(t,(int)ind(j,t));
        for (t = 1; t < hor;  t++)
        {
            for (int t0 = 0; t0 <= t; t0++)
            {
                //compute price
                price(t,t0) = (alpha(j,t) + (prod(j,t0) + sum(stor(j,Range(t0,t-1)))
                    + cons(j,t0)*coef(t0))* beta(j,t)) / (2 * beta(j,t));
                //compute demand
                demand(t,t0) = alpha(j,t) - beta(j,t) * price(t,t0);
                //compute cost
                sumc = 0;
                for (int i = t0; i <= t; i++)
                {
                    sumc += (prod(j,t0) + cons(j,t0)*coef(t0) + sum(stor(j,Range(t0,i-1))) - price(i,t0))*demand(i,t0);
                }
                c(t0) = sumc + setup(j,t0);
            }
            //find minimal criterium
            f(t+1) = min(c(Range(0,t)) + f(Range(0,t)));
            //result for period t and product j
            ind(j,t) = min(minIndex(c(Range(0,t))+ f(Range(0,t))));
            results(j,t)=price(t,(int)ind(j,t));
        }

    ////OUTPUT
        if (verbose >2)
        {
            for (t = 0; t < hor; t++)
            {
                printf("%d\t%d\t%f\t%d\n",j+1,t+1,-f(t+1),(int)ind(j,t)+1);
            }
        }
    ////OUTPUT

    }

////OUPOUT
    if (verbose>2)
    {
        printf("-------------------------------------\n");
    }
////OUTPUT

    return ind;
}




