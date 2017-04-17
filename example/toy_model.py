# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 12:15:06 2015

@author: noore

Description:
    Generate a toy reaction network model with N intermediate metabolites
    and N+1 reactions.

"""

from ecm.cost_function import EnzymeCostFunction
import numpy as np

N = 2

S = np.matrix(np.zeros((N + 2, N + 1)))
for i in xrange(N + 1):
    S[i, i] = -1.0
    S[i+1, i] = 1.0
v        = np.matrix(np.ones((N + 1, 1)))
kcat     = np.matrix(np.ones((N + 1, 1))) * 100
dG0      = -5.0 * np.matrix(np.ones((N + 1, 1)))
K_M      = np.matrix(np.ones(S.shape))
K_M[S < 0] = 3e-2
K_M[S > 0] = 1e-2

lnC_bounds = np.matrix(np.ones((S.shape[0], 2))) * np.log(1e-4)
lnC_bounds[1:-1, 0] = np.log(1e-6)
lnC_bounds[1:-1, 1] = np.log(1e-2)

ecf1 = EnzymeCostFunction(S, v, kcat, dG0, K_M, lnC_bounds,
                          ecf_version='ECF1')
mdf1_sol, params = ecf1.MDF()
print mdf1_sol
lnC0 = params['ln concentrations']
print ecf1.ECM(lnC0).T

ecf2 = EnzymeCostFunction(S, v, kcat, dG0, K_M, lnC_bounds,
                          ecf_version='ECF2')
mdf2_sol, params = ecf2.MDF()
print mdf2_sol
lnC0 = params['ln concentrations']
print ecf2.ECM(lnC0).T

