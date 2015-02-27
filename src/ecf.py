# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 15:40:11 2015

@author: noore
"""

import itertools
import sys
import types
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as font_manager
from scipy.optimize import minimize
from scipy.integrate import ode
from optimized_bottleneck_driving_force import Pathway

RT = 8.31e-3 * 298.15

class ThermodynamicallyInfeasibleError(Exception):

    def __init__(self, y=None):
        if y is None:
            Exception.__init__(self, 'this reaction system has no feasible solutions')
        else:
            Exception.__init__(self, 'y = %s : is thermodynamically infeasible' 
                                     % str(np.exp(y)))
    
    pass

class NonSteadyStateSolutionError(Exception):
    pass

class ECF(object):
    
    def __init__(self, S, v, kcat, dG0, K_M, A_act, A_inh, K_act, K_inh,
                 V=1.0, c_range=(1e-6, 1e-2), c_fixed=1e-4, ecf_version='ECF4'):
        """
            Construct a toy model with N intermediate metabolites (and N+1 reactions)
            
            Arguments:
                S       - stoichiometric matrix
                v       - steady-state fluxes [umol/min]
                kcat    - specific activities [umol/min/mg]
                dG0     - standard Gibbs free energies of reaction [kJ/mol]
                K_M     - Michaelis-Menten coefficients [M]
                
                V       - cell volume [um^3]
                c_range - allowed range for internal metabolite concentration [M]
                c_fixed - concentration of external metabolites [M]
        """
        self.S = S
        self.v = v
        self.kcat = kcat
        self.dG0 = dG0
        
        assert v.shape == (S.shape[1], 1)
        assert kcat.shape == (S.shape[1], 1)
        assert dG0.shape == (S.shape[1], 1)
        assert S.shape == K_M.shape
        assert S.shape == A_act.shape
        assert S.shape == K_act.shape
        assert S.shape == A_inh.shape
        assert S.shape == K_inh.shape
        
        self.S_subs = abs(self.S)
        self.S_prod = abs(self.S)
        self.S_subs[self.S > 0] = 0
        self.S_prod[self.S < 0] = 0
        self.V = V
        self.y_range = map(np.log, c_range)
        self.y_fixed = np.log(c_fixed)
        
        self.subs_denom = np.matrix(np.diag(self.S_subs.T * np.log(K_M))).T
        self.prod_denom = np.matrix(np.diag(self.S_prod.T * np.log(K_M))).T

        # allosteric regulation term
        self.A_act = A_act
        self.A_inh = A_inh
        self.act_denom = np.matrix(np.diag(self.A_act.T * np.log(K_act))).T
        self.inh_denom = np.matrix(np.diag(self.A_inh.T * np.log(K_inh))).T

        self.Nr = self.S.shape[1]
        self.Nc = self.S.shape[0]
        self.Nint = self.Nc - 2
        
        # check that v is a possible steady-state flux solution
        if (np.abs(self.S[1:-1,:] * v) > 1e-6).any():
            raise NonSteadyStateSolutionError('the provided flux "v" is not a '
            'steady-state solution of the provided stoichiometric model')
        
        if ecf_version == 'ECF4':
            self.ECF = self.ECF4
        elif ecf_version == 'ECF3':
            self.ECF = self.ECF3
        elif ecf_version == 'ECF2':
            self.ECF = self.ECF2
        elif ecf_version == 'ECF1':
            self.ECF = self.ECF1
        else:
            raise ValueError('The enzyme cost function %s is unknown' % ecf_version)
        
    def y_to_lnC(self, y):
        """
            The difference between lnC and y, is that lnC includes also the
            concentrations of the fixed (external) metabolites
        """
        if y.shape == (self.Nint, ):
            y = np.matrix(y).T
        assert y.shape[0] == self.Nint

        return np.vstack([np.ones((1, y.shape[1])) * self.y_fixed,
                          y,
                          np.ones((1, y.shape[1])) * self.y_fixed])

    def get_vmax(self, E):
        """
            calculate the maximal rate of each reaction, kcat is in umol/min/mg and 
            E is in gr, so we multiply by 1000
            
            Returns:
                Vmax  - in units of [umol/min]
        """
        return np.multiply(self.kcat, E) * 1e3 # umol/min
    
    def ECF1(self, lnC):
        """
            C - a matrix of metabolite log-concentrations, where the rows are metabolites in the same order
                as in S, and the columns are indices for different points in the metabolite polytope (i.e.
                conditions).
        """
        # lnC is not used for ECF1, except to determine the size of the result
        # matrix.
        # we multiply by 1e-3 to convert the kcat to umol/min/<gr> instead of <mg>
        return np.tile(np.multiply(self.v, 1e-3/self.kcat), (1, lnC.shape[1]))

    def driving_force(self, lnC):
        # calculate the driving force for every reaction in every condition
        assert lnC.shape[0] == self.Nc
        return -np.tile(self.dG0 / RT, (1, lnC.shape[1])) - self.S.T * lnC

    def eta_thermodynamic(self, lnC):
        df = self.driving_force(lnC)
        
        # replace infeasbile reactions with a positive driving force to avoid negative cost in ECF2
        eta_thermo = 1.0 - np.exp(-df)
        
        # set the value of eta to a negative number when the reaction is infeasible
        # so it will be easy to find them, and also calculating 1/x will not return
        # an error
        eta_thermo[df <= 0] = -1.0
        return eta_thermo

    def eta_kinetic(self, lnC):
        kin_subs = np.exp(self.S_subs.T * lnC - np.tile(self.subs_denom, (1, lnC.shape[1])))
        kin_prod = np.exp(self.S_prod.T * lnC - np.tile(self.prod_denom, (1, lnC.shape[1])))
        eta_kin = kin_subs/(1.0 + kin_subs + kin_prod)
        return eta_kin        

    def eta_allosteric(self, lnC):
        kin_act = np.exp(-self.A_act.T * lnC + np.tile(self.act_denom, (1, lnC.shape[1])))
        kin_inh = np.exp(self.A_inh.T * lnC - np.tile(self.inh_denom, (1, lnC.shape[1])))
        eta_kin = 1.0 / (1.0 + kin_act) / (1.0 + kin_inh)
        return eta_kin

    def ECF2(self, lnC):
        """
            lnC - a matrix of metabolite log-concentrations, where the rows are 
                  metabolites in the same order as in S, and the columns are 
                  indices for different points in the metabolite polytope (i.e.
                  conditions).
        """
        ECF2 = np.multiply(self.ECF1(lnC), 1.0/self.eta_thermodynamic(lnC))
    
        # fix the "fake" values that were given in ECF2 to infeasible reactions
        ECF2[ECF2 < 0] = np.nan
        
        return ECF2

    def ECF3(self, lnC):
        """
            lnC - a matrix of metabolite log-concentrations, where the rows 
                  are metabolites in the same order as in S, and the columns
                  are indices for different points in the metabolite polytope
                  (i.e. conditions).
        """
        # calculate the product of all substrates and products for the kinetic term
        return np.multiply(self.ECF2(lnC), 1.0/self.eta_kinetic(lnC))

    def ECF4(self, lnC):
        """
            Add a layer of allosteric activators/inibitors
        """
        return np.multiply(self.ECF3(lnC), 1.0/self.eta_allosteric(lnC))

    def get_fluxes(self, lnC, E):
        v = np.tile(self.get_vmax(E), (1, lnC.shape[1]))
        v = np.multiply(v, self.eta_thermodynamic(lnC))
        v = np.multiply(v, self.eta_kinetic(lnC))
        v = np.multiply(v, self.eta_allosteric(lnC))
        return v
    
    def generate_contour_data(self, n=300):
        rng = np.linspace(self.y_range[0], self.y_range[1], n)
        X0, X1 = np.meshgrid(rng, rng)
        C = np.matrix(list(itertools.product(rng, repeat=2)))
        
        F = self.y_fixed * np.ones((C.shape[0], 1))
        lnC = np.hstack([F, C] + [F]*(self.Nc-3)).T
         
        E = self.ECF(lnC)

        return X0, X1, E
        
    def plot_contour(self, n=300):
        X0, X1, E = self.generate_contour_data(n)

        # calculate the total cost of all enzymes
        Etot = np.sum(E, axis=0)
        Etot = np.array(np.reshape(Etot.T, X0.shape).T)

        fig, axarr = plt.subplots(1, self.Nr+1, figsize=((self.Nr+1)*3,4),
                                  sharey=True)
        axarr[0].set_ylabel(r'$y_1$ [M]')
        for i, ax in enumerate(axarr):
            ax = axarr[i]
            ax.set_xscale('log')
            ax.set_yscale('log')
            if i == self.Nr:
                tmp = Etot
                ax.set_title(r'$\varepsilon$')
            else:
                tmp = np.array(np.reshape(E[i,:].T, X0.shape).T)
                ax.set_title(r'$\varepsilon_%d$' % (i+1))
            ax.contourf(np.exp(X0), np.exp(X1), tmp,
                        locator=ticker.LogLocator(), cmap=plt.cm.PuBu)
            ax.set_xlabel(r'$y_0$ [M]')
        
        y_min = self.ECM()
        ax.plot([y_min[0,0]], [y_min[1,0]], 'xr', label='ECM')
        
        y_mdf = self.MDF()
        ax.plot([y_mdf[0,0]], [y_mdf[1,0]], 'xy', label='MDF')
        ax.legend(loc='upper left')

    def plot_v_vs_y(self, E, n=50, y1=None):
        fig = plt.figure(figsize=(6,6))
        y1 = y1 or self.y_fixed
        y = np.vstack([np.linspace(-6, -2, n)*np.log(10),
                       np.ones((1, n))*y1])

        v = self.get_fluxes(self.y_to_lnC(y), E)
        
        ax = fig.add_subplot(1, 1, 1, xscale='log', yscale='log')
        for i in xrange(self.Nr):
            ax.plot(np.exp(y[0, :]).flat, v[i, :].flat, linewidth=2,
                    alpha=0.5, label=r'$v_%d$' % i)
        ax.legend(loc='best')
        ax.set_xlabel(r'$y_0$ [M]')
        ax.set_ylabel(r'flux')

    def plot_Y_slice(self, n=1000, Yvalue=None):
        fig = plt.figure(figsize=(6,6))
        Yvalue = Yvalue or self.y_fixed

        lnC = np.ones((4, n))*self.y_fixed
        lnC[1,:] = np.linspace(self.y_range[0], self.y_range[1], n)
        lnC[2,:] = Yvalue
        
        #data = np.vstack([self.ECF1(lnC), self.ECF2(lnC), self.ECF3(lnC)])
        data = self.ECF(lnC)
        data = np.vstack([data, data.sum(axis=0)])
        ax = fig.add_subplot(1, 1, 1, xscale='log', yscale='log')
        ax.plot(lnC[1, :], data.T, '-', linewidth=2)
        ax.legend([r'$R_%d$' % (i+1) for i in xrange(3)] + ['Total'], loc='best')
        ax.set_xlabel(r'$y_0$')
        ax.set_ylabel(r'$\varepsilon$')
        ax.set_title('$y_1$ = %.e M' % Yvalue)
        ax.set_ylim(ymax=1e2)
        fig.tight_layout()

    def plot_Y_slices(self, nx=300, ny=6):
        fig = plt.figure(figsize=(8,8))
        x = 10**np.linspace(-4, -2, nx)
        y = 10**np.linspace(-3.6, -2.4, ny)
        X, Y = np.meshgrid(x, y)
        
        C = np.ones((4, nx*ny))*self.c_fixed
        C[1,:] = X.flat
        C[2,:] = Y.flat
        lnC = np.log(C)
        
        E = self.ECF(lnC)
        Etot = np.sum(E, axis=0)
        Etot = np.array(np.reshape(Etot.T, X.shape).T)

        ax = fig.add_subplot(1, 1, 1, xscale='log', yscale='log')
        ax.plot(x, Etot, '-', linewidth=2)
        ax.legend(map(lambda i:'[B] = %.1e M' % i, y), loc='lower right')
        ax.set_xlabel('$y_0$')
        ax.set_ylabel('$\varepsilon$')
        ax.set_ylim(ymax=1e3)
        fig.tight_layout()

    def ECM(self):
        
        def optfun(y):
            e = self.ECF(self.y_to_lnC(y)).sum(axis=0)[0,0]
            #print lnx.T, e
            if np.isnan(e):
                return 1e5
            else:
                return np.log(e)
                
        def constfun(y):
            return self.driving_force(self.y_to_lnC(y))
        
        y0 = self.MDF()
        bounds = self.y_range

        res = minimize(optfun, y0, bounds=[bounds]*self.Nint, method='TNC')
        if not res.success:
            print res
        return np.matrix(res.x).T

    def MDF(self):
        """
            Find an initial point (x0) for the optimization using MDF.
        """
        y_bounds = np.array([(self.y_fixed, self.y_fixed)] + \
                            [self.y_range] * self.Nint + \
                            [(self.y_fixed, self.y_fixed)])
        
        p = Pathway(self.S, self.v, self.dG0/RT, y_bounds[:, 0], y_bounds[:, 1])
        mdf, params = p.FindMDF()
        if np.isnan(mdf) or mdf < 0.0:
            raise ThermodynamicallyInfeasibleError()
        return np.log(params['concentrations'][1:-1, 0])
    
    def is_feasible(self, y):
        lnC = self.y_to_lnC(y)
        df = self.driving_force(lnC)
        return (df > 0).all()

    def simulate(self, E, y0=None, t_max=10, dt=0.1, eps=1e-5, figure=None):
        """
            Find the steady-state solution for the metabolite concentrations
            given the enzyme abundances
            
            Arguments:
                E    - enzyme abundances [gr]
                y0   - initial concentration of internal metabolites (default: MDF solution)
                eps  - the minimal change under which the simulation will stop
            
            Returns:
                v    - the steady state flux
                y    - the steady state internal metabolite concentrations
        """
        if type(E) == types.ListType:
            assert len(E) == self.Nr
            E = np.matrix(E).T
        else:
            assert E.shape == (self.Nr, 1)
            
        if (E <= 0).any():
            return np.nan, np.nan
        
        # normalize the total enzyme amount to 1000 mg
        #E *= (1000.0/E.sum(0))

        def f(t, y):
            # we only care about the time derivatives of the internal metabolites
            # (i.e. the first and last one are assumed to be fixed in time)
            return self.S[1:-1,:] * self.get_fluxes(self.y_to_lnC(y), E)
        
        if y0 is None:
            y0 = np.ones((self.Nint, 1)) * self.y_fixed
        else:
            assert y0.shape == (self.Nint, 1)

        if not self.is_feasible(y0):
            raise ThermodynamicallyInfeasibleError(y0)
        
        v = self.get_fluxes(self.y_to_lnC(y0), E)

        T = np.array([0])
        Y = y0.T
        V = v.T
        
        r = ode(f)
        r.set_initial_value(y0, 0)
        
        while r.successful() and \
              r.t < t_max and \
              (r.t < 0.05*t_max or (np.abs(self.S[1:-1,:] * v) > eps).any()):
            r.integrate(r.t + dt)
            v = self.get_fluxes(self.y_to_lnC(r.y), E)

            T = np.hstack([T, r.t])
            Y = np.vstack([Y, r.y.T])
            V = np.vstack([V, v.T])
        
        if r.t >= t_max:
            v_inf = np.nan
            y_inf = np.nan
        else:
            v_inf = V[-1, 0]
            y_inf = Y[-1,:]

        if figure is not None:
            prop = font_manager.FontProperties(size=10)

            ax1 = figure.add_axes([0.1, 0.2, 0.2, 0.6], frameon=True)
            ax1.set_yscale('log')
            ax1.plot(T, np.exp(Y))
            ax1.set_xlabel('time [min]')
            ax1.set_ylabel('concentration [M]')
            ax1.set_ylim(np.exp(self.y_range))
            ax1.legend(map(lambda i: '$y_%d$'% i, xrange(self.Nint)), loc='best',
                       frameon=False, prop=prop)
            
            X0, X1, E = self.generate_contour_data()
            Etot = np.sum(E, axis=0)
            Etot = np.array(np.reshape(Etot.T, X0.shape).T)
    
            ax2 = figure.add_axes([0.4, 0.2, 0.2, 0.6], frameon=True)
            ax2.plot(T, V)
            ax2.set_xlabel('time [min]')
            ax2.set_ylabel('flux [$\mu$mol/min]')
            ax2.set_yscale('log')
            ax2.legend(map(lambda i: '$v_%d$'% i, xrange(self.Nr)), loc='best',
                       frameon=False, prop=prop)
    
            ax3 = figure.add_axes([0.7, 0.2, 0.2, 0.6], frameon=True)
            ax3.set_xscale('log')
            ax3.set_yscale('log')
            ax3.contourf(np.exp(X0), np.exp(X1), Etot,
                         locator=ticker.LogLocator(),
                         cmap=plt.cm.PuBu)
            ax3.plot(np.exp(Y[0,0]), np.exp(Y[0,1]), 'x', markersize=5,
                     color='r', label=r'$y(0)$')
            ax3.plot(np.exp(Y[:,0]), np.exp(Y[:,1]), '.', markersize=3,
                     color='y', label=r'$y(t)$')
            ax3.set_xlabel('$y_0$ [M]')
            ax3.set_ylabel('$y_1$ [M]')

            if np.isfinite(y_inf).all():
                ytmp = list(np.exp(y_inf).flat)
                ax3.plot(ytmp[0], ytmp[1], 'x', markersize=5,
                         color='g', label=r'$y(\infty)$')
                ax3.set_title(r'$y(\infty) = [%s]$' % 
                              ','.join(map(lambda x : '%.1e' % np.exp(x), y_inf.flat)),
                              fontsize=12)
            ax3.legend(loc='upper left', frameon=False, numpoints=1, prop=prop)

        return v_inf, y_inf
        
    def simulate_3D(self, y0=None, n=30, figure=None):
        if self.Nr != 3:
            raise Exception('You can only use simulate_3D for models of 3 enzymes')
            
        if y0 is None:
            y0 = np.ones((self.Nint, 1)) * self.y_fixed
        if not self.is_feasible(y0):
            raise ThermodynamicallyInfeasibleError('initial point y0 = %s' % str(y0))
        
        if self.Nint != 2:
            raise ValueError('can only simulate ECM for a network with 2 internal metabolites')

        # map points on the unit square to the L2 unit sphere.
        # (0, 0) -> (0, 1, 0)
        # (0, 1) -> (1, 0, 0)
        # (1, 0) -> (0, 0, 1)
        # this has to be an affine transformation (since 0 is not mapped to 0)
        # we add another 'fake' dimension which is always 1
        proj2to3 = np.matrix([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]) * \
                   np.matrix([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]]).I
        
        rng = np.linspace(0.0, 1.0, n+3)
        
#        emat = np.matrix(np.ones((3, n*(n+1)/2)))
#        k = 0
#        for i in xrange(n):
#            for j in xrange(n-i):
#                emat[0, k] = rng[i]
#                emat[1, k] = rng[j]
#                k += 1
        emat = np.matrix(list(itertools.product(rng, repeat=2))).T
        emat = np.vstack([emat, np.ones((1, emat.shape[1]))])
        E = proj2to3 * emat
        
        V = []
        print 'Simulating dynamic system for multiple enzyme concentrations ...'
        for i in xrange(E.shape[1]):
            sys.stderr.write('%d%%\r' % (i * 100.0 / E.shape[1]))
            if (E[:, i] <= 0).any():
                V.append(0)
                continue
            v, _ = self.simulate(np.abs(E[:, i:i+1]), y0, eps=1e-3)
            
            # normalize the flux by the total amount of the enzyme
            # since we are maximizing the flux per enzyme
            V.append(v / float(E[:, i:i+1].sum(0)))
        print '[DONE]'
        V = np.array(V)
        if np.isnan(V).all():
            raise Exception('None of the simulations converged')
        i_max = np.nanargmax(V)
        
        E_max = E[:, i_max:i_max+1]
        v_max, y_max = self.simulate(E_max, y0)
        E_over_v_max = E_max * (self.v[0, 0] / v_max) # rescale E to match the enzyme cost per flux self.v
        
        if figure is not None:
            ax3 = figure.add_subplot(1, 1, 1)
            ax3.set_title(r'$v / \varepsilon$ [umol/min/g]')
            V = V.reshape((n+3, n+3))
            
            proj2to2 = np.matrix([[-1.0, 1.0, 0.0], [0.0, 0.0, 2.0]]) * \
                       np.matrix([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]]).I
            X = proj2to2 * emat
            x = X[0, :].reshape((n+3, n+3))
            y = X[1, :].reshape((n+3, n+3))
            mappable = ax3.contourf(x, y, V,
                                   cmap=plt.cm.PuBu)
            plt.colorbar(mappable, ax=ax3)
            
            ## mark the maximal steady-state flux point
            ax3.plot(X[0, i_max], X[1, i_max], 'y.',
                     label=r'$\varepsilon_{max} = (%s)$' % 
                           ','.join(map(lambda x : '%.2f' % x, E_max.flat)))
            ax3.legend(loc='upper left', numpoints=1, frameon=False)

            ## plot the 3 pure enzyme distributions (orthogonal basis)
            ax3.plot([-1, 1, 0, -1], [0, 0, 2, 0], 'r-')
            ax3.text(-1, -0.1, r'$\varepsilon = (0,1,0)$')
            ax3.text(0.9, -0.1, r'$\varepsilon = (1,0,0)$')
            ax3.text(0.1, 2, r'$\varepsilon = (0,0,1)$')
            ax3.set_xlim(-1.5, 1.5)
            ax3.set_ylim(-0.5, 2.5)
        
        return E_over_v_max, y_max, E
    
    @staticmethod
    def _make_figure(s, E):
        fig = plt.figure(figsize=(12, 4))
        fig.suptitle(r'%s, $\varepsilon$ = <%s>, $\varepsilon_{\rm tot}$ = %.2f [mg]' % \
                     (s,
                     ','.join(map(lambda x : '%.1f' % (x*1e3), E.flat)),
                     1e3*E.sum(0)), fontsize=10)
        return fig
        
    def generate_pdf_report(self, pdf_fname):
        y_mdf = self.MDF()
        E_mdf = self.ECF(self.y_to_lnC(y_mdf))

        y_ecm = self.ECM()
        E_ecm = self.ECF(self.y_to_lnC(y_ecm))

        fig1 = ECF._make_figure('MDF', E_mdf)
        self.simulate(E_mdf, figure=fig1)

        fig2 = ECF._make_figure('ECM', E_ecm)
        self.simulate(E_ecm, figure=fig2)
        
        fig4 = plt.figure(figsize=(6, 5))
        E_max, y_max, E = self.simulate_3D(30, figure=fig4)
        
        fig3 = ECF._make_figure(r'$\max\left(v / \sum{\varepsilon}\right)$', E_max)
        self.simulate(E_max, figure=fig3)

        pp = PdfPages(pdf_fname)
        pp.savefig(fig1)
        pp.savefig(fig2)
        pp.savefig(fig3)
        pp.savefig(fig4)
        pp.close()
        