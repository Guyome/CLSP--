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
    verbose = _verbose;
    alpha = new Array<double,2>(_alpha, shape(product,period), neverDeleteData);
    beta = new Array<double,2>(_beta, shape(product,period), neverDeleteData);
    prod = new Array<double,2>(_prod, shape(product,period), neverDeleteData);
    stor = new Array<double,2>(_stor, shape(product,period), neverDeleteData);
    cons = new Array<double,2>(consumption, shape(product,period), neverDeleteData);
    setupcost = new Array<double,2>(_setup, shape(product,period), neverDeleteData);
    setup = new Array<double,2>(product,period);
    price = new Array<double,2>(product,period);
    production = new Array<double,2>(product,period);
    storage = new Array<double,2>(product,period);
    ind = new Array<int,2>(product,period);
    constraint = new Array<double,1>(_constraint, shape(period), neverDeleteData);
    coef = new Array<double,1>(period);
    //initial point
    (*coef) = 0.;//KKT initiated as null
    (*storage) = 0;
    (*production) = 0;
    (*price) = alpha->copy();
    
 ////OUTPUT
    if (verbose >2)
    {
        HeurClsp::plotParam();
    }
///OUTPUT
}

void HeurClsp::plotVariables()
{
    printf("\nJ\tT\tPrice\t\tProd.\t\tHold.\t\tSetup\t\tCoef\n");
    printf("----------------------------------------------------------------\n");
    for(int j = 0; j < product; j ++)
    {
        for(int t = 0; t < period; t ++)
        {
            printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n",j,t,(*price)(j, t),(*production)(j, t),(*storage)(j, t),(*setup)(j, t),(*coef)(t));
        }
         printf("----------------------------------------------------------------\n");
    }
}

void HeurClsp::plotParam()
{
    printf("\nJ\tT\tSlope\t\tInter.\t\tProd. C.\tHold. C.\tSetup C.\tCons. per P.\tConst.\n");
    printf("----------------------------------------------------------------------------------------------\n");
    for (int j = 0; j < product; j ++)
    {
        for (int t = 0; t < period; t ++)
        {
            printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t\t%f\n",j,t,(*alpha)(j,t),(*beta)(j,t),(*prod)(j,t),(*stor)(j,t),(*setupcost)(j,t),(*cons)(j,t),(*constraint)(t));
        }
        printf("----------------------------------------------------------------------------------------------\n");
    }
}

double** HeurClsp::ArrayToPtr(Array<double,2> array)
{
    double** ptr;
    ptr = new double*[product];
    for (int j = 0; j < product; j ++)
    {
        ptr[j] = new double[period];
        for (int t = 0; t < period; t ++)
        {
            ptr[j][t] = array(j,t);
        }
    }
    return ptr;
}

double** HeurClsp::getPrice()
{
    return ArrayToPtr((*price));
}

double** HeurClsp::getProd()
{
    return ArrayToPtr((*production));
}

double** HeurClsp::getHold()
{
    return ArrayToPtr((*storage));
}

double** HeurClsp::getSetup()
{
    return ArrayToPtr((*setup));
}

double* HeurClsp::getCoef()
{
    double* ptr;
    ptr = new double[period];
    for (int t = 0; t < period; t ++)
    {
            ptr[t] = (*coef)(t);
    }
    return ptr;
}

void HeurClsp::thomas()
{
    Array<double,1> c(period);//cost function
    Array<double,1> f(period+1);//productective function
    Array<double,2> tprice(period,period);//local price 
    Array<double,2> tdemand(period,period);//local demand
    double sumc;
    int t;//current time period

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
        (*ind)(j,t) = 0;
        (*price)(j,t)=tprice(t,(int)(*ind)(j,t));
        (*production)(j,t) = (*alpha)(j,t)-(*beta)(j,t)*(*price)(j,t);
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
            if ( t > (*ind)(j,t) )
            {
                (*storage)(j,(*ind)(j,t)) += (*alpha)(j,t) - (*beta)(j,t)*(*price)(j,t);
            }
        }
        (*setup) = where((*production) > 0,1,0);
    }

////OUPOUT
    if (verbose >2)
    {
        HeurClsp::plotVariables();
    }
////OUTPUT
}

void HeurClsp::coefheur()
{
    Array<double,1> consValue(period);
    Array<double,1> sortlist(product);
    Array<int,1> linkedproduct(period);
    int obj,tps,t0,counttps,countobj,countt0,nblink;

    //find the first violated constraint
    for (int t = 0; t < period; t ++)
    {
        consValue(t) = sum( (*cons)(Range::all(),t)*(*production)(Range::all(),t));
    }
    //initiate to zero
    tps = first( consValue > (*constraint) );
    counttps = 0;
    while ( (tps >= 0) & (counttps < period))
    {

    ////OUPOUT
        if (verbose >2)
        {
            printf("\nConstraint violated at time %d (%f > %f)\n",tps,consValue(tps),(*constraint)(tps));
        }
    ////OUPOUT

        //while the constraint is violated
        //remove all production for the obj product
        //order product per consumption
        sortlist = (*cons)(Range::all(), tps)*(*setup)(Range::all(), tps);
        obj = 0;
        countobj = 0; //initiate counter
        while ( (consValue(tps) > (*constraint)(tps)) & (obj >= 0) & (countobj < product))
        {        
            //selecte the product with the bigest consumption
            obj = max(maxIndex(sortlist));
            //compute the number of time period who are produce in tps
            linkedproduct = (*ind)(obj, Range::all()).copy();
            nblink = count(linkedproduct == tps);
            //search all time period who are producted in tps (see thomas)
            t0 = tps;
            countt0 = 0;
            while ( (consValue(tps) > (*constraint)(tps)) & (countt0 < nblink))
            {
                //no demand for the obj product at period tps
                (*production)(obj, tps) -= (*alpha)(obj, t0) - (*beta)(obj, t0)*(*price)(obj, t0);
                //update storage if we have hold something
                if ( t0 > tps )
                {
                    (*storage)(obj, tps) -= (*alpha)(obj, t0) - (*beta)(obj, t0)*(*price)(obj, t0);
                }
                //price out
                (*price)(obj, t0) = (*alpha)(obj, t0)/(*beta)(obj, t0);
                //update constraint value
                consValue(tps) = sum( (*cons)(Range::all(),tps)*(*production)(Range::all(),tps));

            ////OUPOUT
                if (verbose >2)
                {
                    printf("\tIn period %d: drop the object number %d producted for the period %d\n",tps,obj,t0);
                }
            ////OUTPUT
                //update loop stoping condition
                linkedproduct(t0) = -1;
                t0 = first(linkedproduct == tps);
                countt0 ++;
            }
            //count the number of loop
            sortlist(obj) = -1;
            countobj ++;
        }
        //modify the obj's production to saturate the tps constraint
        (*coef)(tps) = ( (*constraint)(tps) 
            - sum( (*cons)(Range::all(), tps)*(*production)(Range::all(), tps) ) )
            / (*cons)(obj, tps);
        (*production)(obj, tps) = (*coef)(tps);
        (*price)(obj, tps) = ( (*alpha)(obj, tps) - (*production)(obj, tps) ) / (*beta)(obj, tps);
        consValue(tps) = sum( (*cons)(Range::all(),tps)*(*production)(Range::all(),tps));

    ////OUPOUT
        if (verbose >2)
        {
            printf("\tNew constraint value for period %d: %f\n",tps,consValue(tps));
        }
    ////OUPOUT
        
        //find the next violated constraint
        tps = first( consValue > (*constraint) );
        //count the number of loop
        counttps ++;
    }

////OUPOUT
    if (verbose >2)
    {
        HeurClsp::plotVariables();
    }
////OUTPUT

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
            - (*setup)(j,Range::all())*(*setupcost)(j,Range::all())
            );
    }
    return sum(profit);
}

bool HeurClsp::feasible()
{
    Array<double,1> consValue(period);
    for (int t = 0; t < period; t ++)
    {
        consValue(t) = sum( (*cons)(Range::all(),t)*(*production)(Range::all(),t));
    }
    return all( consValue <= (*constraint) );
}

double HeurClsp::heursolver()
{
    Array<double,1> previouscoef(period);
    double diff, upper, lower;
    int count = 0;
    firstIndex t;
    diff = eps + 1.;

    while ( (diff > eps) & (count < cycle) )
    {

    ////OUPOUT
        if (verbose >2)
        {
            printf("\nITER\t diff\n");
            printf("--------------------\n");
            printf("%d\t\t%f\n",count,diff);
        }
    ////OUPOUT

        previouscoef = (*coef);
        //compute lower and upper bound
        thomas();
        upper = objective();
        //add relaxation to the objective
        for (int t = 0; t < period; t ++)
        {
            upper += (*coef)(t)*( (*constraint)(t)-sum( (*cons)(Range::all(),t)*(*production)(Range::all(),t) ) );
        }
        //stop iteration if upper bound is feasible
        if (HeurClsp::feasible())
        {
            lower = upper;
            count = cycle;
        }
        else
        {
            coefheur();
            lower = objective();
            //update KKT coefficients
            (*coef) = param*previouscoef + (1-param)*(*coef);
            //update stoping conditions 
        }
        diff = upper - lower;
        count ++;
    }
////OUPOUT
    if (verbose >1)
    {
        printf("\nLast objective:\t\t\t %f\n",lower);
        printf("Number of iteration:\t\t\t %d\n",count);
        printf("Difference between upper and lower bound:\t %f\n",diff);
    }
////OUTPUT

    return lower;
}


