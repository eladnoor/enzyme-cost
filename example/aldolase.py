# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 11:43:14 2015

@author: noore

Description:
    create a toy diamon network model (i.e. two parallel pathways from A to D, via
    intermediates B or C)

"""

import numpy as np
from ecm.cost_function import EnzymeCostFunction
import matplotlib.pyplot as plt

R = 8.31e-3
DEFAULT_TEMP = 298.15  # K
RT = R * DEFAULT_TEMP

Nr = 3
Nc = 4

S = np.matrix(np.zeros((Nc, Nr)))

S[0, 0] = -1
S[1, 0] = 1
S[2, 0] = 1
S[1, 1] = -1
S[2, 1] = 1
S[2, 2] = -1
S[3, 2] = 1

v          = np.matrix([1.0, 1.0, 2.0]).T
kcat       = np.matrix([1.0, 1.0, 1.0]).T
dGm_r      = np.matrix([-3.0, -2.0, -3.0]).T
dG0_r      = dGm_r - RT * S.T * np.matrix(np.ones((Nc, 1))) * np.log(1e-3)
K_M        = np.matrix(np.ones(S.shape))
K_M[S < 0] = 9e-2
K_M[S > 0] = 1e-2

lnC_bounds = np.log(np.matrix([[1e-9]*Nc, [1e-1]*Nc]).T)

A_act      = np.matrix(np.zeros(S.shape))
A_inh      = np.matrix(np.zeros(S.shape))
K_act      = np.matrix(np.ones(S.shape))
K_inh      = np.matrix(np.ones(S.shape))

toy_ecf = EnzymeCostFunction(S, v, kcat, dG0_r, K_M, lnC_bounds, None, None,
                             A_act, A_inh, K_act, K_inh)


#%%
kcat_list = np.logspace(-5, 5, 100)
costs_list = []
for kcat in kcat_list:
    toy_ecf.kcat[1] = kcat
    lnC = None
    while lnC is None:
        lnC = toy_ecf.ECM()
    costs = toy_ecf.ECF(lnC)
    costs_list.append(list(costs.flat))
        
costs_list = np.array(costs_list)
#%%
fig1 = plt.figure(figsize=(5, 5))
ax = fig1.add_subplot(1, 1, 1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(kcat_list, costs_list)
ax.plot(kcat_list, costs_list.sum(1))
ax.legend(['enzyme 0 cost', 'enzyme 1 cost', 'enzyme 2 cost', 'total cost'])
ax.set_xlabel('$k_{cat}$ of enzyme 1')

