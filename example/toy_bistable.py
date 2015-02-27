# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 15:01:27 2015

@author: noore
"""

from ecf import ECF
import numpy as np

import os
if not os.path.exists('../res'):
    os.mkdir('../res')

N = 2

S = np.matrix(np.zeros((N + 2, N + 1)))
for i in xrange(N + 1):
    S[i, i] = -1.0
    S[i+1, i] = 1.0
v        = np.matrix(np.ones((N + 1, 1)))  # [umol/min]
kcat     = np.matrix([40e-0, 10e-0, 40e-0]).T # [umol/min/mg]
dG0      = np.matrix([-50, 0, -50]).T      # [kJ/mol]
K_M      = np.matrix(np.ones(S.shape))     # [M]
K_M[S < 0] = 1e-9
K_M[S > 0] = 1e-9

# add a positive allosteric feedback from compound 2 to reaction 0
A_act    = np.matrix(np.zeros(S.shape)) 
K_act    = np.matrix(np.ones(S.shape))
A_act[2, 0] = 1
K_act[2, 0] = 1e-4
A_act[1, 2] = 1
K_act[1, 2] = 1e-2

A_inh    = np.matrix(np.zeros(S.shape)) 
K_inh    = np.matrix(np.ones(S.shape)) 

ecf = ECF(S, v, kcat, dG0, K_M, A_act, A_inh, K_act, K_inh)

E = np.matrix([1e-3, 1e-3, 1e-3]).T # [g]

Y = [[3e-6, 2e-6], [3e-5, 2e-5], [3e-4, 2e-4], [3e-3, 2e-3],
     [3e-3, 1e-6], [3e-3, 1e-5], [3e-3, 1e-4], [3e-3, 1e-3]]

for i in xrange(len(Y)):
    y0 = np.log(np.matrix(Y[i])).T
    y_str = ','.join(map(lambda x : '%.1e' % x, Y[i]))
    fig = ECF._make_figure('ECM : $y_0$ = <%s> [M]' % y_str, E)
    ecf.simulate(E, y0=y0, figure=fig, eps=1e-2, dt=0.01, t_max=10)

#ecf.generate_pdf_report('../res/bistable.pdf')

