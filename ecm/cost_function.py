# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 15:40:11 2015

@author: noore
"""

import numpy as np
from scipy.optimize import minimize
from optimized_bottleneck_driving_force import Pathway
from component_contribution.thermodynamic_constants import default_RT as RT
from util import CastToColumnVector

class NonSteadyStateSolutionError(Exception):
    pass

class ThermodynamicallyInfeasibleError(Exception):

    def __init__(self, y=None):
        if y is None:
            Exception.__init__(self, 'this reaction system has no feasible solutions')
        else:
            Exception.__init__(self, 'y = %s : is thermodynamically infeasible' 
                                     % str(np.exp(y)))
    pass

class EnzymeCostFunction(object):

    ECF_LEVEL_NAMES = ['capacity [M]', 'thermodynamic', 'saturation', 'allosteric']
    
    def __init__(self, S, flux, kcat, dG0, KMM, lnC_bounds,
                 A_act=None, A_inh=None, K_act=None, K_inh=None,
                 ecf_version='ECF4'):
        """
            Construct a toy model with N intermediate metabolites (and N+1 reactions)
            
            Arguments:
                S        - stoichiometric matrix [unitless]
                v        - steady-state fluxes [M/s]
                kcat     - turnover numbers [1/s]
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
        
        assert flux.shape       == (self.Nr, 1)
        assert kcat.shape       == (self.Nr, 1)
        assert dG0.shape        == (self.Nr, 1)
        assert KMM.shape        == (self.Nc, self.Nr)
        assert lnC_bounds.shape == (self.Nc, 2)

        self.S = S
        self.flux = flux
        self.kcat = kcat
        self.dG0 = dG0
        self.KMM = KMM
        self.ecf_version = ecf_version
        
        self.S_subs = abs(self.S)
        self.S_prod = abs(self.S)
        self.S_subs[self.S > 0] = 0
        self.S_prod[self.S < 0] = 0

        # convert the ranges to logscale        
        self.lnC_bounds = lnC_bounds
        
        # allosteric regulation term

        if A_act is None or K_act is None:
            self.A_act = np.matrix(np.zeros(S.shape))
            self.K_act = np.matrix(np.ones(S.shape))
        else:
            assert S.shape == A_act.shape
            assert S.shape == K_act.shape
            self.A_act = A_act
            self.K_act = K_act

        if A_inh is None or K_inh is None:
            self.A_inh = np.matrix(np.zeros(S.shape))
            self.K_inh = np.matrix(np.ones(S.shape))
        else:
            assert S.shape == A_inh.shape
            assert S.shape == K_inh.shape
            self.A_inh = A_inh
            self.K_inh = K_inh
        
        # preprocessing: these auxiliary matrices help calculate the ECF3 and 
        # ECF4 faster.
        self.subs_denom = np.matrix(np.diag(self.S_subs.T * np.log(self.KMM))).T
        self.prod_denom = np.matrix(np.diag(self.S_prod.T * np.log(self.KMM))).T
        self.act_denom = np.matrix(np.diag(self.A_act.T * np.log(self.K_act))).T
        self.inh_denom = np.matrix(np.diag(self.A_inh.T * np.log(self.K_inh))).T

        try:        
            self.ECF = eval('self.' + ecf_version)
        except AttributeError:
            raise ValueError('The enzyme cost function %s is unknown' % ecf_version)
        
    def Serialize(self):
        mdict = {'S' : self.S,
                 'flux': self.flux,
                 'kcat': self.kcat,
                 'dG0': self.dG0,
                 'lnC_bounds': self.lnC_bounds,
                 'KMM': self.KMM,
                 'A_act': self.A_act,
                 'A_ihn': self.A_inh,
                 'K_act': self.K_act,
                 'K_inh': self.K_inh,
                 'ecf_version': self.ecf_version}
        return mdict
        
    def _DrivingForce(self, lnC):
        """
            calculate the driving force for every reaction in every condition
        """
        assert lnC.shape[0] == self.Nc
        if len(lnC.shape) == 1:
            return -self.dG0 / RT - np.dot(self.S.T, lnC)
        else:
            return -np.tile(self.dG0 / RT, (1, lnC.shape[1])) - self.S.T * lnC

    def _EtaThermodynamic(self, lnC):
        assert lnC.shape[0] == self.Nc
        df = self._DrivingForce(lnC)
        
        # replace infeasbile reactions with a positive driving force to avoid negative cost in ECF2
        eta_thermo = 1.0 - np.exp(-df)
        
        # set the value of eta to a negative number when the reaction is infeasible
        # so it will be easy to find them, and also calculating 1/x will not return
        # an error
        eta_thermo[df <= 0] = -1.0
        return eta_thermo

    def _STerm(self, lnC):
        assert lnC.shape[0] == self.Nc
        return np.exp(self.S_subs.T * lnC - np.tile(self.subs_denom, (1, lnC.shape[1])))

    def _PTerm(self, lnC):
        assert lnC.shape[0] == self.Nc
        return np.exp(self.S_prod.T * lnC - np.tile(self.prod_denom, (1, lnC.shape[1])))

    def _EtaKinetic(self, lnC):
        """
            the kinetic part of eMC3(1SP), i.e. S / (1 + S + P)
        """
        S = self._STerm(lnC)
        P = self._PTerm(lnC)
        return S/(1.0 + S + P)

    def _EtaAllosteric(self, lnC):
        assert lnC.shape[0] == self.Nc
        kin_act = np.exp(-self.A_act.T * lnC + np.tile(self.act_denom, (1, lnC.shape[1])))
        kin_inh = np.exp(self.A_inh.T * lnC - np.tile(self.inh_denom, (1, lnC.shape[1])))
        eta_kin = 1.0 / (1.0 + kin_act) / (1.0 + kin_inh)
        return eta_kin

    def IsFeasible(self, lnC):
        assert lnC.shape == (self.Nc, 1)
        df = self._DrivingForce(lnC)
        return (df > 0).all()
        
    def GetVmax(self, E):
        """
            calculate the maximal rate of each reaction, kcat is in umol/min/mg and 
            E is in gr, so we multiply by 1000
            
            Returns:
                Vmax  - in units of [umol/min]
        """
        assert E.shape == (self.Nr, 1)
        return np.multiply(self.kcat, E) # in M/s
    
    def ECF1(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector
        
            Returns:
                The most basic Enzyme Cost Function (only dependent on flux
                and kcat). Gives the predicted enzyme concentrations in [M]
        """
        # lnC is not used for ECF1, except to determine the size of the result
        # matrix.
        assert lnC.shape == (self.Nc, 1)
        return np.tile(np.multiply(self.flux, 1.0/self.kcat), (1, lnC.shape[1]))

    def ECF2_S(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector
        
            Returns:
                The thermodynamic-only Enzyme Cost Function.
                Gives the predicted enzyme concentrations in [M].
        """
        assert lnC.shape == (self.Nc, 1)
        ECF2 = np.multiply(self.ECF1(lnC), 1.0/self._EtaThermodynamic(lnC))
    
        # fix the "fake" values that were given in ECF2 to infeasible reactions
        ECF2[ECF2 < 0] = np.nan
        
        return ECF2

    def ECF3_1SP(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector
        
            Returns:
                An Enzyme Cost Function that integrates kinetic and thermodynamic
                data, but no allosteric regulation.
                Gives the predicted enzyme concentrations in [M].
        """
        # calculate the product of all substrates and products for the kinetic term
        assert lnC.shape == (self.Nc, 1)
        return np.multiply(self.ECF2_S(lnC), 1.0/self._EtaKinetic(lnC))

    def ECF4(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector
        
            Returns:
                The full Enzyme Cost Function, i.e. with kinetic, thermodynamic
                and allosteric data.
                Gives the predicted enzyme concentrations in [M].
        """
        assert lnC.shape == (self.Nc, 1)
        return np.multiply(self.ECF3SP(lnC), 1.0/self._EtaAllosteric(lnC))

    def GetEnzymeCostPartitions(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector
        
            Returns:
                A matrix contining the enzyme costs separated to the 4 ECF
                factors (as columns).
                The first column is the ECF1 predicted concentrations in [M].
                The other columns are unitless (added cost, always > 1)
        """
        assert lnC.shape == (self.Nc, 1)
        cap = self.ECF1(lnC)                  # capacity
        trm = 1.0/self._EtaThermodynamic(lnC) # thermodynamics
        kin = 1.0/self._EtaKinetic(lnC)       # kinetics
        alo = 1.0/self._EtaAllosteric(lnC)    # allostery
        return np.hstack([cap, trm, kin, alo])

    def GetFluxes(self, lnC, E):
        assert len(lnC.shape) == 2
        assert lnC.shape[0] == self.Nc
        assert E.shape == (self.Nr, 1)
        
        v = np.tile(self.GetVmax(E), (1, lnC.shape[1]))
        v = np.multiply(v, self._EtaThermodynamic(lnC))
        v = np.multiply(v, self._EtaKinetic(lnC))
        v = np.multiply(v, self._EtaAllosteric(lnC))
        return v
    
    def ECM(self, lnC0=None):
        """
            Use convex optimization to find the y with the minimal total
            enzyme cost per flux, i.e. sum(ECF(lnC))
        """
        
        def optfun(lnC):
            lnC = CastToColumnVector(lnC)
            minimal_df = self._DrivingForce(lnC).min()
            if minimal_df > 0: # all reactions are feasible
                e = self.ECF(lnC).sum(axis=0)[0,0]
                if np.isnan(e):
                    raise Exception('ECF returns NaN although all reactions are feasible')
                return np.log(e)
            else: # give a large penalty proportional to the negative driving force
                return 1e20 * abs(minimal_df)
                
        assert lnC0.shape == (self.Nc, 1)

        bounds = zip(self.lnC_bounds[:,0].flat, self.lnC_bounds[:,1].flat)
        res = minimize(optfun, lnC0, bounds=bounds, method='SLSQP')
        #res = minimize(optfun, lnC0, bounds=bounds, method='TNC')
        
        if not res.success:
            print res
        
        lnC_min = np.matrix(res.x).T
        return lnC_min

    def MDF(self):
        """
            Find an initial point (x0) for the optimization using MDF.
        """
        p = Pathway(self.S, self.flux, self.dG0, self.lnC_bounds)
        return p.FindMDF()
    
