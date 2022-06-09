#!/usr/bin/python
"""
@author: Boris Deshev
This will calculate bins with (semi-)equal number of points

It requires three parameters to be provided by the user:
1.Table name
2.Name of the column which will be binned
3.Number of bins

24.Jul.2018
"""

import astropy
from astropy.table import *#Table
import numpy as np
import sys
np.random.seed(0)

def check_arguments(arguments):
    #Check if all three parameters are provided
    while (len(arguments) != 3) and (arguments[-1] != 'q') and (arguments[-1] != 'Q'):
        arguments = input("Please give Table name, Column name, number of bins\nOr pres q to exit!\n").split()
    if (arguments[-1] == 'q') or (arguments[-1] == 'Q'):
        print('***********\nI hope you are not disappointed!')
        sys.exit()
    calc_bins(arguments)
    
def calc_bins(arguments):
    #Read in the input table
    if type(arguments[0]) is not astropy.table.table.Table:
        try:
            tbl = Table.read(arguments[0])
        except ValueError:
            format = raw_input('Cannot read the table. What format is it?\n')
            tbl = Table.read(arguments[0], format=format)
    
    else:
        tbl = arguments[0]
    tn = float(len(tbl))
    #Check if the table is longer than the required number of bins
    if tn <= float(arguments[2]):
        print('**********\nThe table has fewer entries than the required number of bins!')
        print('Table length is', tn)
        print('The number of bins is', arguments[2])
        sys.exit()

    tbl.sort(keys=arguments[1])
    minval = np.min(tbl[arguments[1]])
    maxval = np.max(tbl[arguments[1]])
    #Check if the column contains a range of values and no NaNs
    if minval == maxval:
        print('The minimum and maximum are the same!')
        sys.exit()

    #The number of points per bin will be
    binpop = int(tn/float(arguments[2]))
    print('Each bin will contain approximately', binpop, 'points')
    
    """
    Fill in the list with the borders of the bins
    It may happen that the left over from the rounding above adds up to
    a very nig number of points in the last bin. This will try to scatter them randomly
    without significantly affecting the population number of indivdual bins
    """
    indices = np.arange(0, len(tbl), int(len(tbl)/arguments[2]))[:arguments[2]]

    leftover = tn%arguments[2]
    if leftover != 0:
        add_to_bins = np.random.choice(range(arguments[2]), int(leftover)-1, replace=False)

        for toadd in add_to_bins:
            indices[toadd:] += 1

    bins = tbl[arguments[1]][indices].data

    # If the uper bound of the last bin is not higher than the max data value 
    # np.digitize would just invent anothe bin (crappy isn't it)
    bins = np.append(bins, maxval+1)
    # Guess what, it does exactly the same with the low boundary too
    bins[0] -= 1
    
    #Formulate the output
    print('*'*20)
    print('RESULTS:')
    print('The bin edges are:')
    print(["%.3e" %ff for ff in bins])
    
    print('\nThe number of points in each bin is:')
    print(np.histogram(tbl[arguments[1]],bins=bins)[0])
    
    return bins

if __name__ == "__main__":
    check_arguments(sys.argv)

