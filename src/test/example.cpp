#include "../HeurClsp.hpp"
#include <stdio.h>

int main(int argv, char* argc[])
{
    int verb = 2;
    int J = 1;
    int T = 4;
    int cycle = 100;
    double param = 0.3;
    double eps = 1.;
    double a[] = {100.,100.,100.,100.};//slope demand function
    double b[] = {1.,1.,1.,1.};//intercept demand function
    double v[] = {20.,20.,20.,20.};//production cost
    double h[] = {2.,2.,2.,2.};//storage cost
    double r[] = {18, 20., 32., 12.};//constraint
    double c[] = {1.,1.,1.,1.};//product consumption
    double s[] = {1,1,1,1};//setup
    
    HeurClsp* example = new HeurClsp(a,b,v,h,c,s,r,T,J,verb,cycle,eps,param);

    printf("OK\n");
}
