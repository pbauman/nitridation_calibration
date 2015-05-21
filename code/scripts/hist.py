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
from scipy.stats import gaussian_kde,skew,skewtest,kurtosis
import h5py

h5_file = h5py.File('chain_data.h5', 'r')

print "Reading Data"
gamma_CN_data = h5_file['/constant_gamma_n_constant_gamma_cn/run1/Metropolis-Hastings/log10_gamma_CN_raw_chain']
gamma_N_data = h5_file['/constant_gamma_n_constant_gamma_cn/run1/Metropolis-Hastings/log10_gamma_N_raw_chain']

gamma_CN_sub_data = gamma_CN_data[100000:4000000:20]
gamma_N_sub_data = gamma_N_data[100000:4000000:20]

print "Computing stats"
gamma_CN_mean = mean(gamma_CN_sub_data)
gamma_CN_sigma = std(gamma_CN_sub_data)

gamma_N_mean = mean(gamma_N_sub_data)
gamma_N_sigma = std(gamma_N_sub_data)

print "Gamma N stats:"
print "Mean: ", gamma_N_mean
print "Std. Dev.: ", gamma_N_sigma
print "Skewness: ", skew(gamma_N_sub_data)
print "Kurtosis: ", kurtosis(gamma_N_sub_data)
print " "

print "Gamma CN stats:"
print "Mean: ", gamma_CN_mean
print "Std. Dev.: ", gamma_CN_sigma
print "Skewness: ", skew(gamma_CN_sub_data)
print "Kurtosis: ", kurtosis(gamma_CN_sub_data)
print " "

def gaussian(mu,sigma,x):
    return 1.0/(sigma*sqrt(2*pi))*exp(-0.5*((x-mu)/sigma)**2);


tick_label_fontsize=18
axis_label_fontsize=18
matplotlib.rc('xtick', labelsize=tick_label_fontsize )
matplotlib.rc(('xtick.major','xtick.minor'),  pad=10)
matplotlib.rc('ytick', labelsize=tick_label_fontsize)

fig = plot.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel( r"$\log_{10}\gamma_{CN}$", fontsize=axis_label_fontsize)
ax1.set_ylabel( "Output Samples", fontsize=axis_label_fontsize)

print "Plotting Histogram"
ax1.hist( gamma_CN_sub_data, bins=250, align='mid' )

ax2 = ax1.twinx()

#print "Computing KDE"
#density = gaussian_kde(gamma_CN_data)

x = linspace(-3.25, -3.05, 500)

print "Plotting KDE"
#ax2.plot( x, density(x), 'k--', linewidth=1 )
ax2.plot(x, gaussian(gamma_CN_mean,gamma_CN_sigma,x), 'r-', linewidth=1 )

ax2.set_xbound(lower=-3.15, upper=-3.05)
ax2.set_ylabel( "PDF", fontsize=axis_label_fontsize)
#ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
#ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0) )

ax1.grid(True)

print "Saving Figure"
plot.savefig("gamma_CN.pdf", bbox_inches='tight')

############################################
# Average N
############################################

fig2 = plot.figure()
ax3 = fig2.add_subplot(111)
ax3.set_xlabel( r"$\log_{10}\gamma_{N}$", fontsize=axis_label_fontsize)
ax3.set_ylabel( "Output Samples", fontsize=axis_label_fontsize)

print "Plotting Histogram"
ax3.hist( gamma_N_sub_data, bins=250, align='mid' )

ax4 = ax3.twinx()

#print "Computing KDE"
#density = gaussian_kde(gamma_CN_data)

x = linspace(-4.5, -4.1, 500)

print "Plotting KDE"
#ax2.plot( x, density(x), 'k--', linewidth=1 )
ax4.plot(x, gaussian(gamma_N_mean,gamma_N_sigma,x), 'r-', linewidth=1 )

ax4.set_xbound(lower=-4.5, upper=-4.1)
ax4.set_ylabel( "PDF", fontsize=axis_label_fontsize)
#ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
#ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0) )

ax3.grid(True)

print "Saving Figure"
plot.savefig("gamma_N.pdf", bbox_inches='tight')

# plot.figure(2)
# plot.xlabel( r"$\gamma_{N}$", fontsize=axis_label_fontsize)
# plot.ylabel( "Output Samples", fontsize=axis_label_fontsize)

# p3 = plot.hist( data[:,1], bins=250, align='mid' )

# ax2 = plot.twinx()

# density = gaussian_kde(data[:,1])

# #x2 = linspace(0.4,0.6,500)

# #p4 = plot.plot( x2, density(x2), 'k-', linewidth=2 )

# plot.ylabel( "PDF", fontsize=axis_label_fontsize)
# ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0) )

# plot.figure(4)
# p5 = plot.acorr(data[:,1])

# plot.figure(3)
# plot.xlabel( r"$T_a$", fontsize=axis_label_fontsize)
# plot.ylabel( "Output Samples", fontsize=axis_label_fontsize)

# p5 = plot.hist( data[:,2], bins=50, align='mid' )

# ax3 = plot.twinx()

# density = gaussian_kde(data[:,2])

# x3 = linspace(0.5, 2.0,100)

# p6 = plot.plot( x3, density(x3), 'k-', linewidth=2 )

# plot.ylabel( "PDF", fontsize=axis_label_fontsize)
# ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0) )


print "Showing figure"
plot.show()
