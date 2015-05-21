#!/usr/bin/python

import sys
import os
import glob

import matplotlib
from matplotlib import rc
rc('text',usetex=True)

# Different modes
# By default, you can "show()" the figure which gives you an interactive window
# and it will save a .png when you call savefig().
# PDF produces a .pdf file, but show() doesn't seem to work.
#matplotlib.use("PDF")

import matplotlib.pyplot as plot
import numpy as np
from scipy.stats import gaussian_kde,skew,skewtest,kurtosis
import h5py

h5_file = h5py.File('chain_data.h5', 'r')

print "Reading Data"
gamma_CN_data = h5_file['/constant_gamma_n_constant_gamma_cn/run1/Metropolis-Hastings/log10_gamma_CN_raw_chain']
gamma_N_data = h5_file['/constant_gamma_n_constant_gamma_cn/run1/Metropolis-Hastings/log10_gamma_N_raw_chain']

gamma_CN_sub_data = gamma_CN_data[100000:4000000:20]
gamma_N_sub_data = gamma_N_data[100000:4000000:20]

data = np.zeros((gamma_CN_sub_data.size, 2))
data[:,0] = gamma_CN_sub_data
data[:,1] = gamma_N_sub_data

np.savetxt("datapoints.dat",data)
