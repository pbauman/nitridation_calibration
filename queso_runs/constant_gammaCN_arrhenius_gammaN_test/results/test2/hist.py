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
fig1, ax1 = plot.subplots()
ax1.set_xlabel( r"$\gamma_{CN}$", fontsize=axis_label_fontsize)
ax1.set_ylabel( "Output Samples", fontsize=axis_label_fontsize)

ax1.hist( data[:,0], bins=50, align='mid' )

ax2 = ax1.twinx()

density = gaussian_kde(data[:,0])

x = linspace(2, 6,100)

p2 = ax2.plot( x, density(x), 'k-', linewidth=2 )

ax2.set_ylabel( "PDF", fontsize=axis_label_fontsize)
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0) ) 
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ) 

#plot.savefig("ablation_rate_pdf_2010.pdf")

fig2, ax3 = plot.subplots()
ax3.set_xlabel( r"$\gamma_{N}$", fontsize=axis_label_fontsize)
ax3.set_ylabel( "Output Samples", fontsize=axis_label_fontsize)

ax3.hist( data[:,1], bins=50, align='mid' )

ax4 = ax3.twinx()

density1 = gaussian_kde(data[:,1])

x2 = linspace(0,2.0e-3,500)

ax4.plot( x2, density1(x2), 'k-', linewidth=2 )

ax4.set_ylabel( "PDF", fontsize=axis_label_fontsize)
ax4.ticklabel_format(style='sci', axis='x', scilimits=(0,0) ) 
ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ) 


fig3, ax5 = plot.subplots()
ax5.set_xlabel( r"$T_a$", fontsize=axis_label_fontsize)
ax5.set_ylabel( "Output Samples", fontsize=axis_label_fontsize)

ax5.hist( data[:,2], bins=50, align='mid' )

ax6 = ax5.twinx()

density2 = gaussian_kde(data[:,2])

x3 = linspace(-10, -6,100)

p6 = ax6.plot( x3, density2(x3), 'k-', linewidth=2 )

ax6.set_ylabel( "PDF", fontsize=axis_label_fontsize)
ax6.ticklabel_format(style='sci', axis='x', scilimits=(0,0) ) 
ax6.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ) 


plot.show()
