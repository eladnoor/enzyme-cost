# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 15:01:27 2015

@author: noore
"""
import os
from ecm.cost_function import EnzymeCostFunction
from ecm.simulator import EnzymeCostSimulator
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

if not os.path.exists('res'):
    os.mkdir('res')

Nc = 3
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

# v0: X0 -> X1
S[0, 0] = -1.0
S[1, 0] = 1.0
v[0] = 1.0
kcat[0] = 300.0
dG0[0] = -8
K_M[0, 0] = 1.2e-1
K_M[1, 0] = 1.2e-1

# v1: X1 -> X2
S[1, 1] = -2.0
S[2, 1] = 1.0
v[1] = 2.0
kcat[1] = 1000.0
dG0[1] = -100
K_M[1, 1] = 1.2e-2
K_M[2, 1] = 1.2e-2

# v2: X0 -> X1
S[0, 2] = -1.0
S[1, 2] = 1.0
v[2] = 1.0
kcat[2] = 0.01
dG0[2] = -8
K_M[0, 2] = 1.2e-9
K_M[1, 2] = 1.2e-1

# add a positive allosteric feedback from X1 to reaction 0
A_act[1, 0] = 6
K_act[1, 0] = 1.8e-4

v = np.matrix(v).T
kcat = np.matrix(kcat).T
dG0 = np.matrix(dG0).T

lnC_bounds = np.matrix([[1e-4] * Nc, [1e-4] * Nc]).T
lnC_bounds[1, 0] = 1e-6
lnC_bounds[1, 1] = 1e-2
ecf = EnzymeCostFunction(S, v, kcat, dG0, K_M, lnC_bounds, A_act, A_inh, K_act, K_inh)

######################################################################

lnC = np.log(np.tile(np.matrix([1e-4]*Nc).T, (1, 100)))
lnC[1, :] = np.linspace(-6, -2, 100) * np.log(10)

E0 = np.matrix([0.95, 0.045, 0.005]).T # [g]
E1 = np.matrix([0.60, 0.25, 0.15]).T   # [g]

E = np.hstack([E0, E1])

figs = []
for i in range(E.shape[1]):
    v = ecf.GetFluxes(lnC, E[:,i:i+1])
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1, 2, 1, xscale='log', yscale='log')
    plt.plot(np.exp(lnC[1, :]).flat, (v[0, :] + v[2, :]).flat, label=r'$v_0 + v_2$')
    plt.plot(np.exp(lnC[1, :]).flat, v[1, :].flat, label=r'$v_1$')
    ax.set_xlabel(r'$y_1$ [M]')
    ax.set_ylabel('flux [M/s]')
    ax.legend()

    ax = fig.add_subplot(1, 2, 2, xscale='log', yscale='linear')
    plt.plot(np.exp(lnC[1, :]).flat, (ecf.S[1,:] * v).flat)
    ax.set_xlabel(r'$y_1$ [M]')
    ax.set_ylabel(r'$\frac{dy_1}{dt}$ [M/s]')
    figs.append(fig)

simu = EnzymeCostSimulator(ecf)

lnC0 = np.log(np.matrix([1e-4]*Nc)).T

alphas = np.linspace(0, 1, 50)
lnC_steady = []
v_steady = []
for alpha in alphas:
    v_inf, lnC_inf = simu.Simulate(lnC0, E0*(1-alpha) + E1*alpha)
    v_steady.append(v_inf)
    lnC_steady.append(lnC_inf[0, 1])

fig = plt.figure(figsize=(6,8))
ax = fig.add_subplot(2, 1, 1, xscale='linear', yscale='log')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$y_1$ steady state')
ax.plot(alphas, map(np.exp, lnC_steady), '-r')
ax = fig.add_subplot(2, 1, 2, xscale='linear', yscale='log')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$v$ steady state')
ax.plot(alphas, v_steady, '-g')
figs.append(fig)

pp = PdfPages('res/bistable_feedback.pdf')
map(pp.savefig, figs)
pp.close()
