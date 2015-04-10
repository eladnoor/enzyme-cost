# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 15:40:11 2015

@author: noore
"""

import logging
import numpy as np
from scipy.optimize import minimize
from ecm.optimized_bottleneck_driving_force import Pathway
from component_contribution.thermodynamic_constants import default_RT as RT

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

class EnzymeCostFunction(object):
    
    def __init__(self, S, flux, kcat, dG0, KMM,
                 A_act=None, A_inh=None, K_act=None, K_inh=None,
                 c_ranges=None, ecf_version='ECF4'):
        """
            Construct a toy model with N intermediate metabolites (and N+1 reactions)
            
            Arguments:
                S        - stoichiometric matrix [unitless]
                v        - steady-state fluxes [umol/min]
                kcat     - specific activities [umol/min/mg]
                dG0      - standard Gibbs free energies of reaction [kJ/mol]
                KMM      - Michaelis-Menten coefficients [M]
                A_act    - Hill coefficient matrix of allosteric activators [unitless]
                A_inh    - Hill coefficient matrix of allosteric inhibitors [unitless]
                K_act    - affinity coefficient matrix of allosteric activators [M]
                K_inh    - affinity coefficient matrix of allosteric inhibitors [M]
                c_ranges - list of pairs indicating the min and max
                           allowed metabolite concentration [M]
        """
        self.Nc, self.Nr = S.shape
        
        assert flux.shape     == (self.Nr, 1)
        assert kcat.shape     == (self.Nr, 1)
        assert dG0.shape      == (self.Nr, 1)
        assert S.shape        == KMM.shape
        assert c_ranges.shape == (self.Nc, 2)

        self.S = S
        self.flux = flux
        self.kcat = kcat
        self.dG0 = dG0
        
        self.S_subs = abs(self.S)
        self.S_prod = abs(self.S)
        self.S_subs[self.S > 0] = 0
        self.S_prod[self.S < 0] = 0

        # convert the ranges to logscale        
        self.y_bounds = np.log(c_ranges)
        
        self.subs_denom = np.matrix(np.diag(self.S_subs.T * np.log(KMM))).T
        self.prod_denom = np.matrix(np.diag(self.S_prod.T * np.log(KMM))).T

        # allosteric regulation term
        self.A_act = A_act
        self.A_inh = A_inh

        if A_act is None or K_act is None:
            self.act_denom = np.matrix(np.zeros((self.Nr, 1)))
            self.A_act = np.matrix(np.zeros(S.shape))
        else:
            assert S.shape == A_act.shape
            assert S.shape == K_act.shape
            self.act_denom = np.matrix(np.diag(self.A_act.T * np.log(K_act))).T
            self.A_act = A_act

        if A_inh is None or K_inh is None:
            self.inh_denom = np.matrix(np.zeros((self.Nr, 1)))
            self.A_inh = np.matrix(np.zeros(S.shape))
        else:
            assert S.shape == A_inh.shape
            assert S.shape == K_inh.shape
            self.inh_denom = np.matrix(np.diag(self.A_inh.T * np.log(K_inh))).T
            self.A_inh = A_inh

        try:        
            self.ECF = eval('self.' + ecf_version)
        except AttributeError:
            raise ValueError('The enzyme cost function %s is unknown' % ecf_version)
        
    def y_to_lnC(self, y):
        """
            The difference between lnC and y, is that lnC includes also the
            concentrations of the fixed (external) metabolites
        """
        if y.shape == (self.Nc, ):
            y = np.matrix(y).T
        assert y.shape[0] == self.Nc

        return np.vstack([np.ones((1, y.shape[1])) * self.y_fixed,
                          y,
                          np.ones((1, y.shape[1])) * self.y_fixed])

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

    def is_feasible(self, y):
        lnC = self.y_to_lnC(y)
        df = self.driving_force(lnC)
        return (df > 0).all()
        
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
        return np.tile(np.multiply(self.flux, 1e-3/self.kcat), (1, lnC.shape[1]))

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
    
    def ECM(self, y0=None):
        """
            Use convex optimization to find the y with the minimal total
            enzyme cost per flux, i.e. sum(ECF(y))
        """
        
        def optfun(y):
            e = self.ECF(self.y_to_lnC(y)).sum(axis=0)[0,0]
            #print lnx.T, e
            if np.isnan(e):
                return 1e5
            else:
                return np.log(e)
                
        def constfun(y):
            return self.driving_force(self.y_to_lnC(y))
        
        if y0 is None:
            y0 = self.MDF()

        res = minimize(optfun, y0, bounds=self.y_bounds, method='TNC')
        if not res.success:
            print res
        return np.matrix(res.x).T

    def MDF(self):
        """
            Find an initial point (x0) for the optimization using MDF.
        """
        p = Pathway(self.S, self.flux, self.dG0/RT,
                    self.y_bounds[:, 0], self.y_bounds[:, 1])
        mdf, params = p.FindMDF()
        if np.isnan(mdf) or mdf < 0.0:
            logging.error('Negative MDF value: %.1f' % mdf)
            raise ThermodynamicallyInfeasibleError()
        return np.log(params['concentrations'])
    


        