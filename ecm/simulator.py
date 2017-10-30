# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:07:17 2015

@author: noore
"""
import numpy as np
from scipy.integrate import ode
from ecm.errors import ThermodynamicallyInfeasibleError

class EnzymeCostSimulator(object):
    
    def __init__(self, ecf):
        self.ecf = ecf
    
    def Simulate(self, lnC0, E, t_max=1000, dt=1, eps=1e-9):
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
        assert lnC0.shape == (self.ecf.Nc, 1)
        assert E.shape == (self.ecf.Nr, 1)
        
        lnC_bounds_diff = self.ecf.lnC_bounds[:, 1] - self.ecf.lnC_bounds[:, 0]
        idx_fixed = list(np.where(lnC_bounds_diff <= 1e-5)[0].flat)
        idx_non_fixed = list(np.where(lnC_bounds_diff > 1e-5)[0].flat)
            
        def f(t, y):
            # we only care about the time derivatives of the internal metabolites
            # (i.e. the first and last one are assumed to be fixed in time)
            lnC = np.log(np.matrix(y).T)
            v = self.ecf.GetFluxes(lnC, E)
            dy = self.ecf.S * v
            dy[idx_fixed, :] = 0
            return dy
            
        if not self.ecf.IsFeasible(lnC0):
            raise ThermodynamicallyInfeasibleError(lnC0)
        
        v = self.ecf.GetFluxes(lnC0, E)

        T = np.array([0])
        Y = np.exp(lnC0).T
        V = v.T
        
        r = ode(f)
        r.set_initial_value(Y.T, 0)
        
        while r.successful() and \
              r.t < t_max and \
              (r.t < 0.05*t_max or (np.abs(self.ecf.S[idx_non_fixed,:] * v) > eps).any()):
            r.integrate(r.t + dt)
            v = self.ecf.GetFluxes(np.log(r.y), E)

            T = np.hstack([T, r.t])
            Y = np.vstack([Y, r.y.T])
            V = np.vstack([V, v.T])
        
        if r.t >= t_max:
            v_inf = np.nan
            lnC_inf = np.nan
        else:
            v_inf = V[-1, 0]
            lnC_inf = np.log(Y[-1,:])

        return v_inf, lnC_inf
        
