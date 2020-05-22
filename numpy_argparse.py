"""
Assignment for Effective Computing
Jade Sauvé - github: jade-sauve
use numpy and argparse
"""

# Import Packages
import numpy as np 
import argparse as ap 
import os
import sys
import pickle
import argparse

from functools import reduce

# function to find the factors of an integer
def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

# run on fjord or local machine
#path_current = !pwd # only works in ipython
path_current = os.path.abspath('.') 

if path_current.split('/')[1] == 'Users':   
	runpath = '/Users/jadesauve/Documents/Python/scripts/data_analysis_tech'
	toolpath = '/Users/jadesauve/Documents/Python/scripts/python_tools/'
elif path_current.split('/')[1] == 'data1':
	runpath = '/data1/effcom/jsauve/data_analysis_tech'
	toolpath = '/data1/effcom/jsauve/python_tools/'
elif path_current.split('/')[1] == 'home':
	runpath = '/data1/effcom/jsauve/data_analysis_tech'
	toolpath = '/data1/effcom/jsauve/python_tools/'
else:
	print('we have an issue with toolpath')

pht = os.path.abspath(toolpath)
if pht not in sys.path:
    sys.path.append(pht)
from toolbox_cmdline import *

# use argparse
# array sizes, n can't be a prime number
# m = 10
# n = 12

# create the parser object
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--m_int', default=10, type=int)
parser.add_argument('-n', '--n_int', default=12, type=int)

# get the arguments
args = parser.parse_args()

m = args.m_int
n = args.n_int
print('For the numpy functions to work, n needs to be 12')

# select the right factors to make an array with more than 1D
fac = factors(n)
fac = list(fac)
fac.sort()
l = int(len(fac)/2)
s1 = fac[l-1]
s2 = fac[l]

# make some arrays
arr1 = np.array(np.arange(0,m))
arr2 = np.array(np.arange(0,n).reshape(s1,s2))

# look up what methods are available
dir(arr1)

# np.all() and np.any()
if np.all(arr1 > 5):
	print('All of array 1 is bigger than 5')
elif np.any(arr1 > 5):
	print('At least one element of array is bigger than 5. ')
else:
	print('No elements of array 1 is bigger than 5.')

# np.base(), returns original object 
a = arr1[1:4]
a.base is arr1

# np.choose
# select one item per column at the specified row
selection = np.choose([0,1,2,1],arr2)

# np.clip(), restrain entry values between 2 and 7
arr1_clipped = np.zeros(arr1.shape)
np.clip(arr1,2,7,arr1_clipped)

# np.cumprod, cumulative product alng an axis (0 will go down a column)
arr2_cp = np.cumprod(arr2, axis=0)

# np.flatten(), 1D copy
arr2_flat = arr2.flatten()

# create directory
# if we are in the correct folder
if path_current == runpath:
	make_dir('numpy_argparse_output')

# write the array to a pickle file
out_fn = path_current+'/numpy_argparse_output/' + 'pickled_array.pkl'
pickle.dump(arr2, open(out_fn, 'wb')) # 'wb' is for write binary

# read the pickle file 
arr3 = pickle.load(open(out_fn, 'rb')) # 'rb is for read binary




