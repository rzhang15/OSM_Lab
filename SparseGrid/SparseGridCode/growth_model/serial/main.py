#======================================================================
#
#     This routine solves an infinite horizon growth model
#     with dynamic programming and sparse grids
#
#     The model is described in Scheidegger & Bilionis (2017)
#     https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400
#
#     external libraries needed:
#     - IPOPT (https://projects.coin-or.org/Ipopt)
#     - PYIPOPT (https://github.com/xuy/pyipopt)
#     - TASMANIAN (http://tasmanian.ornl.gov/)
#
#     Simon Scheidegger, 11/16 ; 07/17
#======================================================================

import nonlinear_solver_initial as solver     #solves opt. problems for terminal VF
import nonlinear_solver_iterate as solviter   #solves opt. problems during VFI
from parameters import *                      #parameters of model
import interpolation as interpol              #interface to sparse grid library/terminal VF
import interpolation_iter as interpol_iter    #interface to sparse grid library/iteration
import postprocessing as post                 #computes the L2 and Linfinity error of the model

import TasmanianSG                            #sparse grid library
import numpy as np


#======================================================================
# Start with Value Function Iteration

# terminal value function
valnew=[TasmanianSG.TasmanianSparseGrid()]*5

if (numstart==0):
    for itheta in range(ntheta):
        valnew[itheta]=interpol.sparse_grid(n_agents, iDepth, 1)
        valnew[itheta].write("valnew_" + str(theta_range[itheta]) + '_'+str(numstart)+".txt") #write file to disk for restart

# value function during iteration
else:
    for itheta in range(ntheta):
        valnew[itheta].read("valnew_" + str(theta_range[itheta]) + '_'+str(numstart)+".txt") #write file to disk for restart

valold=[TasmanianSG.TasmanianSparseGrid()]*5
valold=valnew

for i in range(numstart, numits):
    valnew=[TasmanianSG.TasmanianSparseGrid()]*5
    for itheta in range(ntheta):
        valnew[itheta]=interpol_iter.sparse_grid_iter(n_agents, iDepth, valold, theta_range[itheta])
        valold[itheta]=valnew[itheta]

    for itheta in range(ntheta):
        valnew[itheta].write("valnew_" + str(theta_range[itheta]) + '_'+str(i+1)+".txt")

#======================================================================
print "==============================================================="
print " "
print " Computation of a growth model of dimension ", n_agents ," finished after ", numits, " steps"
print " "
print "==============================================================="
#======================================================================

# compute errors
avg_err=post.ls_error(n_agents, numstart, numits, No_samples)

# plot graphs
graphs=post.plot_routine(n_agents, valnew, n_agents-1, 50)

#======================================================================
print "==============================================================="
print " "
print " Errors are computed -- see errors.txt"
print " "
print "==============================================================="
#======================================================================
