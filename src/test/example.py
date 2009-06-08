#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import heurclsp as sl
from numpy import array

T = 2
J = 2
alpha = [[100,100],[100,100]]
beta = [[1,1],[1,1]]
prod = [[20,20],[20,20]]
stor = [[2,2],[2,2]]
setup = [[5,5],[5,5]]
cons = [[1,1],[1,1]]
constraint = [18,20]

cycle = 100;
eps = 1.;
param = 0.5;
verbose = 3;

test = sl.heurclsp(alpha, beta, prod,
    stor, cons, setup, constraint, T, J, verbose, cycle, eps, param)
test.heursolver()
