#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import random as rd
import csv
from sys import exit
from math import fabs

def generateparam(nb_obj,time_hor,d_mean=100., d_sigma=10,
    d_ratio=100, c_mean=20, c_sigma=2, h_ratio=4, s_ratio=10,
    cs_mean=1, cs_sigma=.1,r_mean=30,r_sigma=10):
    """
    This function generates parameters for PCLSP solver.
    All paramters are considerated as gaussian.
    Slope and intercept are correlated ( see d_ratio)
    All costs are correlated (see h_ratio and s_ratio)
    """
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

def importdata(address):
    """
    This function import data form cvs file define as 
    in the following example:
    
    "T";4
    "J";2
    "Slope";100;100;100;100
    ;100;100;100;100
    "Intercept";1;1;1;1
    ;1;1;1;1
    "Setup cost";3;3;3;3
    ;3;3;3;3
    "Storage cost";2;2;2;2
    ;2;2;2;2
    "Production cost";20;20;20;20
    ;2;2;2;2
    "Consuption";0.3;0.3;0.3;0.3
    ;2;2;2;2
    "Constraint";35;40;42;24
    """
    try:
        inputcsv = csv.reader(open(address, "r"), delimiter=";", lineterminator="\n", quoting=csv.QUOTE_NONNUMERIC)
    except IOError:
        print "File not exists or is unreadable, please check it."
        exit(1)

    data = tuple() # all data
    item = list() # each tabular
    count = 0
    subcount = 0
    try:
        for row in inputcsv:
            if count < 2 : # read Time period and number of product
                data += (int(row[1]),)
            else :
                item.append(row[1:])
                subcount +=1 
                if subcount == data[1]:
                    data += (item,)
                    item = list()
                    subcount = 0
            count += 1
        if (data[1] > 1):
            data += (item[0],) # manage the last tabular
    except:
        print "File is not well formated, please correct it."
        exit('file %s, line %d' % (address, inputcsv.line_num))
    return data

def exportdata(address,slope, intercept, prodcost, holdcost,
    setupcost, consumption, constraint, nb_obj, time_hor):
    """
    This function import data form cvs file define as 
    in the following example:
    
    "T";4
    "J";2
    "Slope";100;100;100;100
    ;100;100;100;100
    "Intercept";1;1;1;1
    ;1;1;1;1
    "Setup cost";3;3;3;3
    ;3;3;3;3
    "Storage cost";2;2;2;2
    ;2;2;2;2
    "Production cost";20;20;20;20
    ;2;2;2;2
    "Consuption";0.3;0.3;0.3;0.3
    ;2;2;2;2
    "Constraint";35;40;42;24
    """
    try:
        outputcsv = csv.writer(open(address, "w"), delimiter=";", lineterminator="\n", quoting=csv.QUOTE_NONNUMERIC)
    except IOError:
        print "Impossible to write or create the file, please check it."
        exit(1)
    
    outputcsv.writerow(["T",time_hor])
    outputcsv.writerow(["J",nb_obj])
    data = {'Slope' : slope, 'Intercept' : intercept,
        'Setup cost' :setupcost, 'Storage cost' : holdcost,
        'Production cost' : prodcost, 'Consuption' : consumption}
    names = ['Slope', 'Intercept',
        'Setup cost', 'Storage cost',
        'Production cost', 'Consuption']
    for item in names:
        for j in xrange(nb_obj):
            if j == 0:
                row = [item]
            else:
                row = ['']
            row.extend(data[item][j])
            outputcsv.writerow(row)
    row = ["Constraint"]
    row.extend(constraint)
    outputcsv.writerow(row)
    
def outGAP(file_name,pclsp,save=True,file_format="svg"):
    import matplotlib.pyplot as plt
    import numpy as np
    plt.plot(100*np.array(pclsp.gap),"r--")
    plt.title(str(int(pclsp.period))+" period(s) and "+str(int(pclsp.product))+" product(s)")
    plt.ylabel("GAP (%)")
    plt.xlabel("Number of cycle")
    if save:
        plt.savefig(file_name+".svg",format=file_format,transparent=True)
    else:
        plt.show()
