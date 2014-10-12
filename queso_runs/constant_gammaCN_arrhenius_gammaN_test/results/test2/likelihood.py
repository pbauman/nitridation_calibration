#!/usr/bin/python

import sys
import os
import glob

from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import rc


rc('text',usetex=True)

# Different modes
# By default, you can "show()" the figure which gives you an interactive window
# and it will save a .png when you call savefig().
# PDF produces a .pdf file, but show() doesn't seem to work.
#matplotlib.use("PDF")

import matplotlib.pyplot as plot
from numpy import * 
from scipy.stats import gaussian_kde

# Column 0: gamma_CN
# Column 1: gamma_N
# Column 2: T_a
data = loadtxt("raw_chain.dat")

likelihood_data = loadtxt("raw_likelihood.dat")

tick_label_fontsize=18
axis_label_fontsize=18
matplotlib.rc('xtick', labelsize=tick_label_fontsize )
matplotlib.rc(('xtick.major','xtick.minor'),  pad=10)
matplotlib.rc('ytick', labelsize=tick_label_fontsize)


#================================================================
# Plot 2010 Data
#================================================================
fig1, ax1 = plot.subplots()
ax1.set_xlabel( r"$\gamma_{N}$", fontsize=axis_label_fontsize)
ax1.set_ylabel( r"$T_a$", fontsize=axis_label_fontsize)

colors = likelihood_data

s1 = ax1.scatter( data[:,1], data[:,2], c=colors )
plot.colorbar(s1)


fig2, ax2 = plot.subplots()
ax2.set_xlabel( r"$\gamma_{CN}$", fontsize=axis_label_fontsize)
ax2.set_ylabel( r"$T_a$", fontsize=axis_label_fontsize)

colors = likelihood_data

s2 = ax2.scatter( data[:,0], data[:,2], c=colors )
plot.colorbar(s2)


fig3, ax3 = plot.subplots()
ax3.set_xlabel( r"$\gamma_{CN}$", fontsize=axis_label_fontsize)
ax3.set_ylabel( r"$\gamma_{N}$", fontsize=axis_label_fontsize)

colors = likelihood_data

s3 = ax3.scatter( data[:,0], data[:,1], c=colors )
plot.colorbar(s3)


fig4 = plot.figure()
ax4 = fig4.add_subplot(111,projection="3d")
s4 = ax4.scatter( data[:,0], data[:,1], data[:,2], c=colors )
plot.colorbar(s4)

plot.show()
