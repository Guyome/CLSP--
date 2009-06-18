#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import heurclsp as sl
import tools as tl


time_hor = 2
nb_obj = 2 
cycle = 100
eps = 1.
param = 0.8
verbose = 3
address = 'test.csv'

#uncomment following comment to
#import data from csv file 
#"""
slope, intercept, prodcost,\
holdcost, setupcost, consumption,\
constraint  = tl.generateparam (nb_obj,time_hor)

tl.exportdata(address,slope, intercept, prodcost, holdcost,
    setupcost, consumption, constraint, nb_obj, time_hor)
#"""
#tl.importdata(address)

test = sl.heurclsp(slope, intercept, prodcost,
    holdcost, consumption, setupcost, constraint, time_hor, nb_obj, verbose, cycle, eps, param)
test.heursolver()
