#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import heurclsp as sl
import tools as tl

#print information about heurclsp lib
print sl.__doc__

time_hor = 2
nb_obj = 2 
cycle = 100
eps = 1.
param = 0.8
verbose = 3
address = 'test.csv'

#uncomment following comment to
#generate random paramter
"""
slope, intercept, prodcost,\
holdcost, setupcost, consumption,\
constraint = tl.generateparam (nb_obj,time_hor)
"""
#uncomment following comment to
#generate random paramter
"""
tl.exportdata(address,slope, intercept, prodcost, holdcost,
    setupcost, consumption, constraint, nb_obj, time_hor)
"""

#import data from csv file 
time_hor, nb_obj, slope, intercept,\
setupcost, holdcost, prodcost,\
consumption, constraint = tl.importdata(address)

print "Initiate\n"
test = sl.pclsp(slope, intercept, prodcost,
    holdcost, consumption, setupcost, constraint, time_hor, nb_obj, verbose, cycle, eps, param)

#print info and mains functions
print test.PCLSPSolvHeur.__doc__
print test.DynProgSolv.__doc__
print test.noQP.__doc__

print "Solve PCSLP\n"
#test.noQP() #uncomment to use heuristic 
test.PCLSPSolvHeur()

#uncomment to print optimals variables
"""
print "Price:","\n\t".join([str(test.price[j]) for j in range(nb_obj)])
print "Production:","\n\t".join([str(test.prod[j]) for j in range(nb_obj)])
print "Holding:","\n\t".join([str(test.hold[j]) for j in range(nb_obj)])
print "Setup:","\n\t".join([str(test.setup[j]) for j in range(nb_obj)])
print "KKT coef:",str(test.kkt)
"""
