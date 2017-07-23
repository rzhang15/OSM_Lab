#======================================================================
#
#     This routine interfaces with the TASMANIAN Sparse grid
#     The crucial part is
#
#     aVals[iI]=solver.initial(aPoints[iI], n_agents)[0]
#     => at every gridpoint, we solve an optimization problem
#
#     Simon Scheidegger, 11/16 ; 07/17
#======================================================================

import TasmanianSG
import numpy as np
from parameters import *
import nonlinear_solver_initial as solver

#======================================================================

def sparse_grid(n_agents, iDepth):

    grid  = TasmanianSG.TasmanianSparseGrid()

    k_range=np.array([k_bar, k_up])

    ranges=np.empty((n_agents, 2))

    for i in range(n_agents):
        ranges[i]=k_range

    iDim=n_agents

    grid1.makeLocalPolynomialGrid(iDim, iOut, iDepth, which_basis, "localp")
    grid1.setDomainTransform(ranges)

    aPoints=grid1.getPoints()
    iNumP1=aPoints.shape[0]
    aVals=np.empty([iNumP1, 1])

    file=open("comparison0.txt", 'w')

    for iK in range(refinement_level):
        grid1.setSurplusRefinement(fTol, 1, "fds")   #also use fds, or other rules
        aPoints = grid1.getNeededPoints()
        aVals = np.empty([aPoints.shape[0], 1])
        for iI in range(iNumP1):
            aVals[iI]=solver.initial(aPoints[iI], n_agents)[0]
            v=aVals[iI]*np.ones((1,1))
            to_print=np.hstack((aPoints[iI].reshape(1,n_agents), v))
            np.savetxt(file, to_print, fmt='%2.16f')

    file.close()
    grid1.loadNeededPoints(aVals)

    aRes = grid1.evaluateBatch(aPnts)
    fError1 = max(np.fabs(aRes[:,0] - aTres))

    print(" {0:9d} {1:9d}  {2:1.2e}".format(iK+1, grid1.getNumPoints(), fError1))

    # write coordinates of grid to a text file
    f2=open("Adaptive_sparse_grid.txt", 'a')
    np.savetxt(f2, aPoints, fmt='% 2.16f')
    f2.close()

    grid2 = TasmanianSG.TasmanianSparseGrid()
    grid2.makeLocalPolynomialGrid(iDim, iOut, refinement_level+iDepth, which_basis, "localp")
    a = grid2.getNumPoints()

    print("\n-------------------------------------------------------------------------------------------------")
    print "   a fix sparse grid of level ", refinement_level+iDepth, " would consist of " ,a, " points"
    print("\n-------------------------------------------------------------------------------------------------\n")

    f=open("grid.txt", 'w')
    np.savetxt(f, aPoints, fmt='% 2.16f')
    f.close()

    return grid
#======================================================================
