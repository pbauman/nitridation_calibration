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
from numpy import * 
from scipy.stats import gaussian_kde

data = loadtxt("raw_chain.dat")

tick_label_fontsize=18
axis_label_fontsize=18
matplotlib.rc('xtick', labelsize=tick_label_fontsize )
matplotlib.rc(('xtick.major','xtick.minor'),  pad=10)
matplotlib.rc('ytick', labelsize=tick_label_fontsize)


#================================================================
# Plot 2010 Data
#================================================================
plot.figure(1)
plot.xlabel( r"$\gamma_{CN}$", fontsize=axis_label_fontsize)
plot.ylabel( "Output Samples", fontsize=axis_label_fontsize)

p1 = plot.hist( data[:,0], bins=50, align='mid' )

ax = plot.twinx()

density = gaussian_kde(data[:,0])

x = linspace(0.6, 0.9,100)

p2 = plot.plot( x, density(x), 'k-', linewidth=2 )

plot.ylabel( "PDF", fontsize=axis_label_fontsize)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0) ) 
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ) 

#plot.savefig("ablation_rate_pdf_2010.pdf")

plot.figure(2)
plot.xlabel( r"$\gamma_{N}$", fontsize=axis_label_fontsize)
plot.ylabel( "Output Samples", fontsize=axis_label_fontsize)

p3 = plot.hist( data[:,1], bins=50, align='mid' )

ax2 = plot.twinx()

density = gaussian_kde(data[:,1])

x2 = linspace(0.1,10,500)

p4 = plot.plot( x2, density(x2), 'k-', linewidth=2 )

plot.ylabel( "PDF", fontsize=axis_label_fontsize)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0) ) 
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ) 


plot.figure(3)
plot.xlabel( r"$T_a$", fontsize=axis_label_fontsize)
plot.ylabel( "Output Samples", fontsize=axis_label_fontsize)

p5 = plot.hist( data[:,2], bins=50, align='mid' )

ax3 = plot.twinx()

density = gaussian_kde(data[:,2])

x3 = linspace(0.2, 1.2,100)

p6 = plot.plot( x3, density(x3), 'k-', linewidth=2 )

plot.ylabel( "PDF", fontsize=axis_label_fontsize)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0) ) 
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ) 


plot.show()
