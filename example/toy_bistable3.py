# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 15:01:27 2015

@author: noore
"""

from ecf import ECF, RT
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import os
if not os.path.exists('../res'):
    os.mkdir('../res')

Nc = 4
Nr = 3

S        = np.matrix(np.zeros((Nc, Nr)))
v        = [0]*Nr                              # [umol/min]
kcat     = [0]*Nr                              # [umol/min/mg]
dG0      = [0]*Nr                              # [kJ/mol]
K_M      = np.matrix(np.ones((Nc, Nr)))        # [M]
A_act    = np.matrix(np.zeros((Nc, Nr)))
K_act    = np.matrix(np.ones((Nc, Nr)))        # [M]
A_inh    = np.matrix(np.zeros((Nc, Nr)))
K_inh    = np.matrix(np.ones((Nc, Nr)))        # [M]

# X0 -> X1 + X2
S[0, 0] = -1.0
S[1, 0] = 0.5
S[2, 0] = 0.5
v[0] = 1.0
kcat[0] = 50.0
dG0[0] = -8
K_M[0, 0] = 1e-1
K_M[1, 0] = 1e-1
K_M[2, 0] = 1e-1

# X1 + X2 -> X3
S[1, 1] = -0.5
S[2, 1] = -0.5
S[3, 1] = 1.0
v[1] = 2.0
kcat[1] = 100.0
dG0[1] = -8
K_M[1, 1] = 1e-1
K_M[2, 1] = 1e-1
K_M[3, 1] = 1e-1

# X0 -> X1 + X2
S[0, 2] = -1.0
S[1, 2] = 0.5
S[2, 2] = 0.5
v[2] = 1.0
kcat[2] = 2e-2
dG0[2] = -8
K_M[0, 2] = 1e-9
K_M[1, 2] = 1e-1
K_M[2, 2] = 1e-1

# add a positive allosteric feedback from X1 to reaction 0
A_act[1, 0] = 4
K_act[1, 0] = 1e-4

v = np.matrix(v).T
kcat = np.matrix(kcat).T
dG0 = np.matrix(dG0).T

ecf = ECF(S, v, kcat, dG0, K_M, A_act, A_inh, K_act, K_inh)

######################################################################
pp = PdfPages('../res/bistable.pdf')

E = np.matrix([0.82, 0.14, 0.05]).T # [g]

lnC = np.log(np.tile(np.matrix([1e-4, 1, 1, 1e-4]).T, (1, 100)))
lnC[1, :] = np.linspace(-6, -2, 100) * np.log(10)
lnC[2, :] = np.linspace(-6, -2, 100) * np.log(10)

v = ecf.get_fluxes(lnC, E)
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1, xscale='log', yscale='log')
plt.plot(np.exp(lnC[1, :]).flat, (v[0, :] + v[2, :]).flat, label=r'$v_0 + v_2$')
plt.plot(np.exp(lnC[1, :]).flat, v[1, :].flat, label=r'$v_1$')
ax.set_xlabel(r'$y_0$ and $y_1$ [M]')
ax.set_ylabel('flux')
ax.legend()

#ecf.plot_v_vs_y(E, y1=np.log(1.1e-4))

Y = np.matrix(np.linspace(-5, -3, 10) * np.log(10))
Y = np.vstack([Y, Y])
#Y = []


for i in xrange(Y.shape[1]):
    y0 = Y[:,i:i+1]
    y_str = ','.join(map(lambda x : '%.1e' % x, y0.flat))
    fig = ECF._make_figure('ECM : $y_0$ = <%s> [M]' % y_str, E)
    ecf.simulate(E, y0=y0, figure=fig, eps=1e-2, dt=0.01, t_max=30)
    pp.savefig(fig)

fig = plt.figure(figsize=(6, 5))
y0_low = np.ones((2, 1)) * np.log(1e-5)
ecf.simulate_3D(y0, 20, figure=fig)
pp.savefig(fig)

fig = plt.figure(figsize=(6, 5))
y0_high = np.ones((2, 1)) * np.log(1e-3)
ecf.simulate_3D(y0, 20, figure=fig)
pp.savefig(fig)

pp.close()
