#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import heurclsp as sl
import random as rd
from math import fabs

def generateparam(nb_obj,time_hor,d_mean=100., d_sigma=10,
    d_ratio=100, c_mean=20, c_sigma=2, h_ratio=4, s_ratio=10,
    cs_mean=1, cs_sigma=.1,r_mean=30,r_sigma=10):
    slope = [[fabs(rd.gauss(d_mean, d_sigma))
        for t in xrange(time_hor)]
        for j in xrange(nb_obj) ]
    intercept = [[fabs(rd.gauss(slope[j][t]/d_ratio, d_sigma/d_ratio**.5))
        for t in xrange(time_hor)]
        for j in xrange(nb_obj) ]
    prodcost = [[fabs(rd.gauss(c_mean, c_sigma))
        for t in xrange(time_hor)]
        for j in xrange(nb_obj) ]
    holdcost = [[fabs(rd.gauss(prodcost[j][t]/h_ratio, c_sigma/h_ratio**.5))
        for t in xrange(time_hor)]
        for j in xrange(nb_obj) ]
    setupcost = [[fabs(rd.gauss(prodcost[j][t]/s_ratio, c_sigma/s_ratio**.5))
        for t in xrange(time_hor)]
        for j in xrange(nb_obj) ]
    consumption = [[fabs(rd.gauss(cs_mean, cs_sigma))
        for t in xrange(time_hor)]
        for j in xrange(nb_obj) ]    
    constraint = [fabs(rd.gauss(r_mean,r_sigma)) for t in xrange(time_hor)]
        
    return slope,intercept, prodcost,\
    holdcost, setupcost, consumption, constraint
    

T = 2
J = 1 
cycle = 100
eps = 1.
param = 0.5
verbose = 3

alpha, beta, prod,\
stor, setup, cons,\
constraint  = generateparam (J,T)

test = sl.heurclsp(alpha, beta, prod,
    stor, cons, setup, constraint, T, J, verbose, cycle, eps, param)
test.heursolver()
