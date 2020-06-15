#!/usr/bin/python

import os,sys,math


stream0 = open(sys.argv[1],'r')
buffer = stream0.readlines()
for line in buffer:
    if line.startswith("DG_LJ = "):
        elems = line.split()
        DG0 = float( elems[2] ) 
        err0 = float( elems[4] )
        break
stream0.close()

stream1 = open(sys.argv[2],'r')
buffer = stream1.readlines()
for line in buffer:
    if line.startswith("DG_LJ = "):
        elems = line.split()
        DG1 = float( elems[2] )
        err1 = float( elems[4] )
        break
stream1.close()

DG = DG1 - DG0
err = math.sqrt( math.pow(err1,2) + math.pow(err0,2) )

print ("DG_LJ_corr = %8.5f +/- %8.5f kcal/mol " % (DG, err))
