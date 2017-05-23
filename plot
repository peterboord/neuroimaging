#! /usr/bin/env python

import sys
import matplotlib.pyplot as pyplot

for filename in sys.argv[1:]:
   with open(filename,'rt') as sf:
         table = []
         for line in sf: table.append( [float(val) for val in line.split()] )
         table = [ row for row in table if len(row) ] ## remove empty rows
         if len(table[0]) == 1 : pyplot.plot( [y[0] for y in table ] )
         for x in xrange(1,len(table[0])): pyplot.plot( [ y[0] for y in table ], [ y[x] for y in table ] )
pyplot.show()
