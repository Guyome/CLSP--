#include <stdio.h>
#include <blitz/array.h>
#include <algorithm>
#include "HeurClsp.hpp"

using namespace blitz;

HeurClsp::HeurClsp(double* _alpha, double* _beta, double* _prod, double* _stor,
        double* consumption, double* _setup, double* _constraint, int _period,
        int _product, int _verbose, int _cycle, double _eps, double _param)
{
    period = _period;
    product = _product;
    cycle = _cycle;
    eps = _eps;
    param = _param;
    alpha = new Array<double,2>(_alpha, shape(product,period), neverDeleteData);
    beta = new Array<double,2>(_beta, shape(product,period), neverDeleteData);
    prod = new Array<double,2>(_prod, shape(product,period), neverDeleteData);
    stor = new Array<double,2>(_stor, shape(product,period), neverDeleteData);
    cons = new Array<double,2>(consumption, shape(product,period), neverDeleteData);
    setup = new Array<double,2>(_setup, shape(product,period), neverDeleteData);
    price = new Array<double,2>(product,period);
    production = new Array<double,2>(product,period);
    storage = new Array<double,2>(product,period);
    ind = new Array<int,2>(product,period);
    constraint = new Array<double,1>(_constraint, shape(period), neverDeleteData);
    coef = new Array<double,1>(period) ; 
    coef = 0;//KKT initiate as null
    
    verbose = _verbose;
}
        
void HeurClsp::thomas()
{
    Array<double,1> c(period);//cost function
    Array<double,1> f(period+1);//productective function
    Array<double,2> tprice(period,period);//local price 
    Array<double,2> tdemand(period,period);//local demand
    double sumc;
    int t;//current time period

////OUTPUT
    if (verbose >2)
    {
        printf("j\tt\tF(t)\t\tInd\n");
        printf("-------------------------------------\n");
    }
////OUTPUT

    for(int j = 0; j < product;  j++)
    {
        //initiate price,demand since there is no storage
        t = 0;
        tprice(t,t) = ((*alpha)(j,t) + ((*prod)(j,t) + (*cons)(j,t)*(*coef)(t))
            * (*beta)(j,t)) / (2 * (*beta)(j,t) );
        tdemand(t,t) = (*alpha)(j,t) - (*beta)(j,t) * tprice(t,t);
        c(t) = tdemand(t,t)*( tprice(t,t) - (*prod)(j,t))-(*setup)(j,t);
        f(t) = 0;
        f(t+1) = -c(t);
        //result for the first period
        (*ind)(j,0) = 0;
        (*price)(j,t)=tprice(t,(int)(*ind)(j,t));
        for (t = 1; t < period;  t++)
        {
            for (int t0 = 0; t0 <= t; t0++)
            {
                //compute price
                tprice(t,t0) = ((*alpha)(j,t) + ((*prod)(j,t0) + sum((*stor)(j,Range(t0,t-1)))
                    + (*cons)(j,t0)*(*coef)(t0))* (*beta)(j,t)) / (2 * (*beta)(j,t));
                //compute demand
                tdemand(t,t0) = (*alpha)(j,t) - (*beta)(j,t) * tprice(t,t0);
                //compute cost
                sumc = 0;
                for (int i = t0; i <= t; i++)
                {
                    sumc += ((*prod)(j,t0) + (*cons)(j,t0)*(*coef)(t0)
                        + sum((*stor)(j,Range(t0,i-1))) - tprice(i,t0))*tdemand(i,t0);
                }
                c(t0) = sumc + (*setup)(j,t0);
            }
            //find minimal criterium
            f(t+1) = min(c(Range(0,t)) + f(Range(0,t)));
            //update all variable for period t
            (*ind)(j,t) = min(minIndex(c(Range(0,t))+ f(Range(0,t))));
            (*price)(j,t) = tprice(t,(int)(*ind)(j,t));
            (*production)(j,t) = 0.; //initiate production
            (*production)(j,(*ind)(j,t)) += (*alpha)(j,t)-(*beta)(j,t)*(*price)(j,t);
            (*storage)(j,t) = max( 0.,(*production)(j,t) - (*alpha)(j,t)-(*beta)(j,t)*(*price)(j,t) );
        }

    ////OUTPUT
        if (verbose >2)
        {
            for (t = 0; t < period; t++)
            {
                printf("%d\t%d\t%f\t%d\n",j+1,t+1,-f(t+1),(int)(*ind)(j,t)+1);
            }
        }
    ////OUTPUT

    }

////OUPOUT
    if (verbose >2)
    {
        printf("-------------------------------------\n");
    }
////OUTPUT
}

void HeurClsp::coefheur()
{
    Array<double,1> consValue(period);
    secondIndex j;
    double saturated;
    int obj;
    int tps;
    
    //find the first violated constaint
    for (int t = 0; t < period; t ++)
    {
        consValue(t) = sum( (*cons)(Range::all(),t)*(*production)(Range::all(),t));
    }
    tps = first( consValue(Range(0,period))>(*constraint)(Range(0,toEnd)) );
    while ( tps < period )
    {
        obj = 0;
        //while the constraint is violated
        //remove all production for the obj product
        while ( consValue(tps) > (*constraint)(tps) )
        {
            //no demand for the obj product at perdior tps
            (*production)(obj, tps) -= (*alpha)(obj, tps)- (*beta)(obj, tps)*(*price)(obj, tps);
            (*price)(obj, tps) = (*alpha)(obj, tps)/(*beta)(obj, tps);
            (*storage)(obj, tps) = (*production)(obj,tps);
            //update constraint value
            consValue(tps) = sum( (*cons)(Range::all(),tps)*(*production)(Range::all(),tps));
            obj ++;
        }
        //modify the obj's production to saturate the tps constraint
        (*production)(obj, tps) -= (*alpha)(obj, tps)- (*beta)(obj, tps)*(*price)(obj, tps);
        saturated = ( (*constraint)(tps) 
            - sum( (*cons)(Range::all(), tps)*(*production)(Range::all(), tps) ) )
            / (*cons)(obj, tps);
        (*price)(obj, tps) = ( (*alpha)(obj, tps) - saturated ) / (*beta)(obj, tps);
        (*production)(obj, tps) += (*alpha)(obj, tps) + (*alpha)(obj, tps)*(*price)(obj, tps);
        (*storage)(obj, tps) = max( 0.,(*production)(obj,tps) - (*alpha)(obj,tps)-(*beta)(obj,tps)*(*price)(obj,tps) );
        //find the next violated constaint
        tps = first( consValue(Range(tps+1,toEnd))>(*constraint)(Range(tps+1,toEnd)) );
        
    ////OUPOUT
        if (verbose >2)
        {
            printf("In period %d: cancellation of the request to the %dth product\n",tps,obj);
        }
    ////OUTPUT
    }
}

double HeurClsp::objective()
{
    Array<double,1> profit(product);
    for (int j = 0; j < product; j ++)
    {
        profit(j) = sum(
            ((*alpha)(j,Range::all())- (*beta)(j,Range::all())*(*price)(j,Range::all()))
            * (*price)(j,Range::all())
            - (*prod)(j,Range::all())*(*production)(j,Range::all())
            - (*stor)(j,Range::all())*(*storage)(j,Range::all())
            //- (*setup)(j,Range::all())*((*production)(j,Range::all) > 0)
            );
    }
    return sum(profit);
}

double HeurClsp::heursolver()
{
    Array<double,1> previouscoef(period);
    double diff, upper,lower;
    int count = 0;
    diff = eps + 1.;
    while ( (diff > eps) & (count < cycle) )
    {
        previouscoef = (*coef);
        //compute lower and upper bound
        thomas();
        upper = objective();
        coefheur();
        lower = objective();
        //update KKT coefficients
        (*coef) = param*previouscoef - (1-param)*(*coef);
        //update stoping conditions 
        diff = upper - lower;
        count ++;
    }
    
////OUPOUT
    if (verbose >1)
    {
        printf("Last objective:\t\t\t %f\n",lower);
        printf("Number of iteration:\t\t\t %d\n",count);
        printf("Difference between upper and lower bound:\t %f",diff);
    }
////OUTPUT

    return lower;
}


