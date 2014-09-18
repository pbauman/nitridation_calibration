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
plot.figure()
plot.xlabel( r"$\gamma_{CN}$", fontsize=axis_label_fontsize)
plot.ylabel( "Output Samples", fontsize=axis_label_fontsize)

p1 = plot.hist( data[:,2], bins=50, align='mid' )

#ax = plot.twinx()

#density = gaussian_kde(data[0,:])

#x = linspace(1.3e-5, 2.6e-5,100)

#p2 = plot.plot( x, density(x), 'k-', linewidth=2 )

#plot.ylabel( "PDF", fontsize=axis_label_fontsize)
#ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0) ) 
#ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ) 

#plot.savefig("ablation_rate_pdf_2010.pdf")


plot.show()
