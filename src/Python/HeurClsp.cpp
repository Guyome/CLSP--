#include <stdio.h>
#include <boost/python.hpp>
#include <blitz/array.h>
#include <algorithm>
#include "HeurClsp.hpp"
#include "QPSolver.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"

using namespace blitz;
using namespace Ipopt;

HeurClsp::HeurClsp(boost::python::list _alpha, boost::python::list _beta, 
    boost::python::list _prod, boost::python::list _stor, boost::python::list consumption,
    boost::python::list _setup, boost::python::list _constraint, int _period,
    int _product, int _verbose, int _cycle, float _eps, float _param)
{
    period = _period;
    product = _product;
    verbose = _verbose;
    cycle = _cycle;
    eps = (double)_eps;
    param = (double)_param;
    alpha = new Array<double,2>(product,period);
    beta = new Array<double,2>(product,period);
    prod = new Array<double,2>(product,period);
    stor = new Array<double,2>(product,period);
    cons = new Array<double,2>(product,period);
    setupcost = new Array<double,2>(product,period);
    setup = new Array<double,2>(product,period);
    price = new Array<double,2>(product,period);
    production = new Array<double,2>(product,period);
    storage = new Array<double,2>(product,period);
    ind = new Array<int,2>(product,period);
    constraint = new Array<double,1>(period);
    coef = new Array<double,1>(period);
    void (HeurClsp::*updatekkt) ();
    this -> updatekkt = &HeurClsp::coefQP;//pointor to function

    //import from python object
    for (int j = 0; j < product; j ++)
    {
        for (int t = 0; t < period; t ++)
        {
            (*alpha)(j,t) = ( extract<double>(_alpha[j][t]) );
            (*beta)(j,t) = ( extract<double>(_beta[j][t]) );
            (*prod)(j,t) = ( extract<double>(_prod[j][t]) );
            (*stor)(j,t) = ( extract<double>(_stor[j][t]) );
            (*setupcost)(j,t) = ( extract<double>(_setup[j][t]) );
            (*cons)(j,t) = ( extract<double>(consumption[j][t]) );
            (*constraint)(t) = ( extract<double>(_constraint[t]) );
        }
    }

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

boost::python::list HeurClsp::ArrayToList(Array<double,2> array)
{
    boost::python::list row;
    boost::python::list col;
    for (int j = 0; j < product; j ++)
    {
        col = boost::python::list();
        for (int t = 0; t < period; t ++)
        {
            col.append( array(j,t) );
        }
        row.append( col );
    }
    
    return row;
}

boost::python::list HeurClsp::getPrice()
{
    return ArrayToList((*price));
}

boost::python::list HeurClsp::getProd()
{
    return ArrayToList((*production));
}

boost::python::list HeurClsp::getHold()
{
    return ArrayToList((*storage));
}

boost::python::list HeurClsp::getSetup()
{
    return ArrayToList((*setup));
}

boost::python::list HeurClsp::getCoef()
{
    boost::python::list col;
    for (int t = 0; t < period; t ++)
    {
            col.append( (*coef)(t) );
    }
    return col;
}

void HeurClsp::setHeur()
{
    updatekkt = &HeurClsp::coefheur;
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

void HeurClsp::coefQP()
{
    //create an instance of QPSolver
    QPSolver* problem = new QPSolver((*alpha),(*beta),(*prod),(*stor),
        (*cons),(*setup),(*constraint),period,product);
    SmartPtr<TNLP> mynlp = problem;
    //create an instance of the IpoptApplication
    SmartPtr<IpoptApplication> app = new IpoptApplication();
    //initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    //set verbosity and derivative test
    app->Options()->SetIntegerValue("print_level", verbose);
    app->Options()->SetStringValue("derivative_test","second-order");
    //run QPSolver
    status = app->Initialize();
    if (status != Solve_Succeeded) 
    {
        //solve
        status = app->OptimizeTNLP(mynlp);
        if ( status == Solve_Succeeded ) 
        {
            //get back variables
            (*price) = problem -> getPrice();
            (*production) = problem -> getProd();
            (*storage) = problem -> getStor();
            (*coef) = problem -> getCoef();
        }
    }

////OUPOUT
    if (verbose >2)
    {
        HeurClsp::plotVariables();
    }
////OUTPUT
}

void HeurClsp::subproblem()
{
    Array<double,1> diff(period);
    Array<double,2> sort(product,period);
    double increment;
    int tps,maxtps,obj;

    sort = cons -> copy();
    for (int t = 0; t < period; t ++)
    {
        (*production)(Range::all(),t) = (*alpha)(Range::all(),t)- (*beta)(Range::all(),t)*(*price)(Range::all(),t);
        diff(t) = sum((*cons)(Range::all(),t)*(*production)(Range::all(),t)) - (*constraint)(t);
        maxtps = t-1;
        while ( (diff(t) > 0) & (maxtps > 0) )
        {
            tps = last( diff(Range(0,maxtps)) < 0 );
            while ( (diff(t) > 0) & (diff(tps) <= 0) )
            {
                obj = max(maxIndex(sort(Range::all(),tps)));
                increment = min(diff(tps)/(*cons)(obj,tps),diff(t)/(*cons)(obj,t));
                (*production)(obj,tps) += increment;
                (*storage)(obj, tps) += increment;
                (*production)(obj,t) -= increment;
                (*storage)(obj, t) = max(0., (*storage)(obj, t) - increment);
                diff(tps) += increment*(*cons)(obj,tps);
                diff(t) -= increment*(*cons)(obj,t);
                sort(obj,tps) = -1;
            }
            maxtps --;
        }
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
            (*this.*updatekkt)();
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

