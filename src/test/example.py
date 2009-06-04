#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import heurclsp as sl
from numpy import array

T = 3
J = 1
alpha = [[100, 100, 100]]
beta = [[1, 1, 1]]
prod = [[20, 20, 20]]
stor = [[2, 2, 2]]
setup = [[5, 5, 5]]
cons = [[1, 1, 1]]
constraint = [ 30, 30, 30]

cycle = 10;
eps = 1.;
param = 0.5;
verbose = 2;

test = sl.heurclsp(alpha, beta, prod,
    stor, cons, setup, constraint, T, J, verbose, cycle, eps, param)
test.heursolver()
