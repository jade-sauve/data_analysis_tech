# -*- coding: utf-8 -*-
# @Author: jadesauve
# @Date:   2020-04-23 12:59:57
# @Last Modified by:   jadesauve
# @Last Modified time: 2020-04-23 13:50:48
"""
read data from 2017-01-0118.ctd
save lists of depth and temp
plot
"""

# imports
import matplotlib.pyplot as plt
import numpy as np


# path
in_dir = '/Users/jadesauve/Coding/data/CTD/' #MODIFY

# define the input filename
in_fn = in_dir + '2017-01-0118.ctd'

#depth,temp = np.genfromtxt(in_fn,skip_header=572,usecols=[1,2])

# define a counter
i=0
# define lists to hold the data
depth=[]
temp=[]
# open the file
with open(in_fn, 'r', errors='ignore') as f:
	for line in f:
		if i >= 570:
			lst = line.split()
			# add the data to the list
			depth.append(float(lst[1]))
			temp.append(float(lst[2]))
		i+=1

# convert the lists to arrays
temp = np.array(temp)
depth = np.array(depth)

# plot
plt.figure()
plt.plot(temp,-depth)
plt.show()
	
