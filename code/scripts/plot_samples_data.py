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

tick_label_fontsize=18
axis_label_fontsize=18
matplotlib.rc('xtick', labelsize=tick_label_fontsize )
matplotlib.rc(('xtick.major','xtick.minor'),  pad=10)
matplotlib.rc('ytick', labelsize=tick_label_fontsize)

chain_data_file = h5py.File('chain_data.h5', 'r')

print "Reading Data"
gamma_CN_data = chain_data_file['/constant_gamma_n_constant_gamma_cn/run1/Metropolis-Hastings/log10_gamma_CN_raw_chain']
gamma_N_data = chain_data_file['/constant_gamma_n_constant_gamma_cn/run1/Metropolis-Hastings/log10_gamma_N_raw_chain']

gamma_CN_sub_data = gamma_CN_data[100000:4000000:20]
gamma_N_sub_data = gamma_N_data[100000:4000000:20]

model_output = np.loadtxt("modelpoints.dat")

fig1 = plot.figure()
ax1 = fig1.add_subplot(111)

ax1.hist( model_output[:,0], bins=250, align='mid')
ax1.plot( (7.6e-7, 7.6e-7), (0, 4500), 'r-')
ax1.errorbar( 7.6e-7, 2000, xerr=2.0e-8, ecolor='red')
ax1.set_xlabel( r"$\Delta m$ [kg]", fontsize=axis_label_fontsize)
ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig1.savefig("A_1_13a_massloss.pdf", bbox_inches='tight')

fig2 = plot.figure()
ax2 = fig2.add_subplot(111)

ax2.hist( model_output[:,1], bins=250, align='mid')
ax2.plot( (1.33e-3, 1.33e-3), (0, 4000), 'r-')
ax2.errorbar( 1.33e-3, 2000, xerr=0.15*1.33e-3, ecolor='red')
ax2.set_xlabel( r"$\overline{\chi_{N}}$", fontsize=axis_label_fontsize)
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig2.savefig("A_1_13a_avgN.pdf", bbox_inches='tight')

fig3 = plot.figure()
ax3 = fig3.add_subplot(111)

ax3.hist( model_output[:,2], bins=250, align='mid')
ax3.plot( (4.1e-7, 4.1e-7), (0, 4000), 'r-')
ax3.errorbar( 4.1e-7, 2000, xerr=2.0e-8, ecolor='red')
ax3.set_xlabel( r"$\Delta m$ [kg]", fontsize=axis_label_fontsize)
ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig3.savefig("O_2_28_massloss.pdf", bbox_inches='tight')

fig4 = plot.figure()
ax4 = fig4.add_subplot(111)

ax4.hist( model_output[:,3], bins=250, align='mid')
ax4.plot( (5.9e-4, 5.9e-4), (0, 4000), 'r-')
ax4.errorbar( 5.9e-4, 2000, xerr=0.15*5.9e-4, ecolor='red')
ax4.set_xlabel( r"$\overline{\chi_{N}}$", fontsize=axis_label_fontsize)
ax4.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig4.savefig("O_2_28_avgN.pdf", bbox_inches='tight')

fig5 = plot.figure()
ax5 = fig5.add_subplot(111)

ax5.hist( model_output[:,4], bins=250, align='mid')
ax5.plot( (6.1e-7, 6.1e-7), (0, 4000), 'r-')
ax5.errorbar( 6.1e-7, 2000, xerr=2.0e-8, ecolor='red')
ax5.set_xlabel( r"$\Delta m$ [kg]", fontsize=axis_label_fontsize)
ax5.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig5.savefig("Q_2_29_massloss.pdf", bbox_inches='tight')

fig6 = plot.figure()
ax6 = fig6.add_subplot(111)

ax6.hist( model_output[:,5], bins=250, align='mid')
ax6.plot( (1.17e-3, 1.17e-3), (0, 4000), 'r-')
ax6.errorbar( 1.17e-3, 2000, xerr=0.15*1.17e-3, ecolor='red')
ax6.set_xlabel( r"$\overline{\chi_{N}}$", fontsize=axis_label_fontsize)
ax6.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig6.savefig("Q_2_29_avgN.pdf", bbox_inches='tight')

fig7 = plot.figure()
ax7 = fig7.add_subplot(111)

ax7.hist( model_output[:,6], bins=250, align='mid')
ax7.plot( (2.04e-6, 2.04e-6), (0, 4000), 'r-')
ax7.errorbar( 2.04e-6, 2000, xerr=2.0e-8, ecolor='red')
ax7.set_xlabel( r"$\Delta m$ [kg]", fontsize=axis_label_fontsize)
ax7.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig7.savefig("U_2_55_massloss.pdf", bbox_inches='tight')

fig8 = plot.figure()
ax8 = fig8.add_subplot(111)

ax8.hist( model_output[:,7], bins=250, align='mid')
ax8.plot( (2.31e-3, 2.31e-3), (0, 4000), 'r-')
ax8.errorbar( 2.31e-3, 2000, xerr=0.15*2.31e-3, ecolor='red')
ax8.set_xlabel( r"$\overline{\chi_{N}}$", fontsize=axis_label_fontsize)
ax8.ticklabel_format(style='sci', axis='x', scilimits=(0,0) )
fig8.savefig("U_2_55_avgN.pdf", bbox_inches='tight')

plot.show()
