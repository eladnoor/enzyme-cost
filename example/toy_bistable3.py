# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 15:01:27 2015

@author: noore
"""
import os
from ecm.cost_function import EnzymeCostFunction
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
dG0[0] = -6
K_M[0, 0] = 1e-1
K_M[1, 0] = 1e-1
K_M[2, 0] = 1e-1

# X1 + X2 -> X3
S[1, 1] = -0.5
S[2, 1] = -0.5
S[3, 1] = 1.0
v[1] = 2.0
kcat[1] = 100.0
dG0[1] = -6
K_M[1, 1] = 1e-1
K_M[2, 1] = 1e-1
K_M[3, 1] = 1e-1

# X0 -> X1 + X2
S[0, 2] = -1.0
S[1, 2] = 0.5
S[2, 2] = 0.5
v[2] = 1.0
kcat[2] = 2e-2
dG0[2] = -6
K_M[0, 2] = 1e-9
K_M[1, 2] = 1e-1
K_M[2, 2] = 1e-1

# add a positive allosteric feedback from X1 to reaction 0
A_act[1, 0] = 4
K_act[1, 0] = 1e-4

v = np.matrix(v).T
kcat = np.matrix(kcat).T
dG0 = np.matrix(dG0).T
lnC_bounds = np.matrix([[1e-6] * 4, [1e-2] * 4]).T

ecf = EnzymeCostFunction(S, v, kcat, dG0, K_M, lnC_bounds, A_act, A_inh, K_act, K_inh)

E = np.matrix([0.82, 0.14, 0.05]).T # [g]
y0_low = np.matrix([np.log(1.2e-5), np.log(1.2e-5)]).T
y0_high = np.matrix([np.log(4.8e-4), np.log(4.8e-4)]).T


#ecf.plot_contour(y0=y0_low)
#ecf.plot_contour(y0=y0_high)
#E_low = ecf.ECF(ecf.y_to_lnC(ecf.ECM(y0_low)))
#E_high = ecf.ECF(ecf.y_to_lnC(ecf.ECM(y0_high)))

#print ecf.simulate(E, y0_low, eps=1e-2, dt=0.01, t_max=30)
#print ecf.simulate(E, y0_high, eps=1e-2, dt=0.01, t_max=30)
#sys.exit(0)

######################################################################


lnC = np.log(np.tile(np.matrix([1e-4, 1, 1, 1e-4]).T, (1, 100)))
lnC[1, :] = np.linspace(-6, -2, 100) * np.log(10)
lnC[2, :] = np.linspace(-6, -2, 100) * np.log(10)

E = np.matrix([[0.85-x, 0.10+x, 0.05] for x in np.linspace(0.0, 0.1, 2)]).T # [g]

figs = []
y0 = np.matrix([np.log(5e-5), np.log(5e-5)]).T
for i in xrange(E.shape[1]):
    v = ecf.GetFluxes(lnC, E[:,i:i+1])
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1, 1, 1, xscale='log', yscale='log')
    plt.plot(np.exp(lnC[1, :]).flat, (v[0, :] + v[2, :]).flat, label=r'$v_0 + v_2$')
    plt.plot(np.exp(lnC[1, :]).flat, v[1, :].flat, label=r'$v_1$')
    ax.set_xlabel(r'$y_0$ and $y_1$ [M]')
    ax.set_ylabel('flux')
    ax.legend()
    figs.append(fig)

    y_str = ','.join(map(lambda x : '%.1e' % np.exp(x), y0.flat))
    fig = EnzymeCostFunction._make_figure('ECM : $y_0$ = <%s> [M]' % y_str, E[:,i:i+1])
    ecf.simulate(E[:,i:i+1], y0=y0, figure=fig, eps=1e-2, dt=0.1, t_max=3000)
    figs.append(fig)


fig = plt.figure(figsize=(6,6))
E0 = np.matrix([0.90, 0.05, 0.05]).T
E1 = np.matrix([0.70, 0.25, 0.05]).T
ecf.simulate_convex(E0, E1, y0=y0, n=100, eps=1e-2, dt=0.1, t_max=3000, figure=fig)
figs.append(fig)


#ecf.plot_v_vs_y(E, y1=np.log(1.1e-4))

#Y = np.matrix(np.linspace(-4.3, -4.2, 10) * np.log(10))
#Y = np.log(np.matrix(np.linspace(5.07e-5, 5.10e-5, 2)))
#Y = np.vstack([Y, Y])
##Y = []
#
#for i in xrange(Y.shape[1]):
#    y0 = Y[:,i:i+1]
#    y_str = ','.join(map(lambda x : '%.1e' % np.exp(x), y0.flat))
#    fig = ECF._make_figure('ECM : $y_0$ = <%s> [M]' % y_str, E)
#    ecf.simulate(E, y0=y0, figure=fig, eps=1e-2, dt=0.01, t_max=30)
#    figs.append(fig)
#

#fig = plt.figure(figsize=(6, 5))
#ecf.simulate_3D(Y[:,0], n=40, eps=1e-2, dt=0.1, t_max=3000, figure=fig)
#figs.append(fig)

#
#fig = plt.figure(figsize=(6, 5))
#ecf.simulate_3D(Y[:,1], n=40, eps=1e-2, dt=0.1, t_max=3000, figure=fig)
#figs.append(fig)

pp = PdfPages('../res/bistable.pdf')
map(pp.savefig, figs)
pp.close()
