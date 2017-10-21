# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 15:40:11 2015

@author: noore
"""

import numpy as np
from scipy.optimize import minimize
from ecm.optimized_bottleneck_driving_force import Pathway
from ecm.util import RT, ECF_DEFAULTS, CastToColumnVector

QUAD_REGULARIZATION_COEFF = 0.2
METABOLITE_WEIGHT_CORRECTION_FACTOR = 1

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

    def __init__(self, S, flux, kcat,
                 dG0, KMM, lnC_bounds,
                 mw_enz=None, mw_met=None,
                 A_act=None, A_inh=None, K_act=None, K_inh=None,
                 params=None):
        """
            Construct a toy model with N intermediate metabolites (and N+1 reactions)

            Arguments:
                S        - stoichiometric matrix [unitless]
                v        - steady-state fluxes [M/s]
                kcat     - turnover numbers [1/s]
                kcat_source - either 'gmean' or 'fwd'
                dG0      - standard Gibbs free energies of reaction [kJ/mol]
                KMM      - Michaelis-Menten coefficients [M]
                A_act    - Hill coefficient matrix of allosteric activators [unitless]
                A_inh    - Hill coefficient matrix of allosteric inhibitors [unitless]
                K_act    - affinity coefficient matrix of allosteric activators [M]
                K_inh    - affinity coefficient matrix of allosteric inhibitors [M]
                c_ranges - list of pairs indicating the min and max
                           allowed metabolite concentration [M]
        """
        self.params = dict(ECF_DEFAULTS)
        self.params.update(params)
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

        self.S_subs = abs(self.S)
        self.S_prod = abs(self.S)
        self.S_subs[self.S > 0] = 0
        self.S_prod[self.S < 0] = 0

        # if the kcat source is 'gmean' we need to recalculate the
        # kcat_fwd using the formula:
        # kcat_fwd = kcat_gmean * sqrt(kEQ * prod_S(KMM) / prod_P(KMM))
        if self.params['kcat_source'] == 'gmean':
            ln_KMM_prod = np.matrix(np.diag(self.S.T * np.log(self.KMM))).T
            ln_ratio = -ln_KMM_prod - self.dG0/RT
            factor = np.sqrt(np.exp(ln_ratio))
            self.kcat = np.multiply(self.kcat, factor)

        # molecular weights of enzymes and metabolites
        if mw_enz is None:
            self.mw_enz = np.zeros((self.Nr, 1))
        else:
            assert mw_enz.shape == ((self.Nr, 1))
            self.mw_enz = mw_enz

        if mw_met is None:
            self.mw_met = np.zeros((self.Nc, 1))
        else:
            assert mw_met.shape == (self.Nc, 1)
            self.mw_met = mw_met

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
        self.D_S_coeff = np.matrix(np.diag(self.S_subs.T * np.log(self.KMM))).T
        self.D_P_coeff = np.matrix(np.diag(self.S_prod.T * np.log(self.KMM))).T
        self.act_denom = np.matrix(np.diag(self.A_act.T * np.log(self.K_act))).T
        self.inh_denom = np.matrix(np.diag(self.A_inh.T * np.log(self.K_inh))).T

        try:
            self.ECF = eval('self._ECF%d' % self.params['version'])
        except AttributeError:
            raise ValueError('The enzyme cost function %d is unknown' %
                             self.params['version'])

        try:
            self.D = eval('self._D_%s' % self.params['denominator'])
        except AttributeError:
            raise ValueError('The denominator function %s is unknown' %
                             self.params['denominator'])

        self.regularization = self.params['regularization']

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
                 'K_inh': self.K_inh}
        mdict.update(self.params)
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

        # replace infeasbile reactions with a positive driving force to avoid
        # negative cost in ECF2
        eta_thermo = 1.0 - np.exp(-df)

        # set the value of eta to a negative number when the reaction is infeasible
        # so it will be easy to find them, and also calculating 1/x will not return
        # an error
        eta_thermo[df <= 0] = -1.0
        return eta_thermo

    def _D_S(self, lnC):
        """
            return a matrix containing the values of D_S
            i.e. prod(s_i / K_i)^n_i

            each row corresponds to a reaction in the model
            each column corresponds to another set of concentrations (assuming
            lnC is a matrix)
        """
        assert lnC.shape[0] == self.Nc
        return np.exp(self.S_subs.T * lnC - np.tile(self.D_S_coeff, (1, lnC.shape[1])))

    def _D_SP(self, lnC):
        """
            return a matrix containing the values of D_SP
            i.e. prod(s_i / K_i)^n_i + prod(p_j / K_j)^n_j

            each row corresponds to a reaction in the model
            each column corresponds to another set of concentrations (assuming
            lnC is a matrix)
        """
        assert lnC.shape[0] == self.Nc
        return np.exp(self.S_subs.T * lnC - np.tile(self.D_S_coeff, (1, lnC.shape[1]))) + \
               np.exp(self.S_prod.T * lnC - np.tile(self.D_P_coeff, (1, lnC.shape[1])))

    def _D_1S(self, lnC):
        """
            return a matrix containing the values of D_1S
            i.e. 1 + prod(s_i / K_i)^n_i

            each row corresponds to a reaction in the model
            each column corresponds to another set of concentrations (assuming
            lnC is a matrix)
        """
        return 1.0 + self._D_S(lnC)

    def _D_1SP(self, lnC):
        """
            return a matrix containing the values of D_1SP
            i.e. 1 + prod(s_i / K_i)^n_i + prod(p_j / K_j)^n_j

            each row corresponds to a reaction in the model
            each column corresponds to another set of concentrations (assuming
            lnC is a matrix)
        """
        return 1.0 + self._D_SP(lnC)

    def _D_CM(self, lnC):
        """
            return a matrix containing the values of D_CM
            i.e. prod(1 + s_i / K_i)^n_i + prod(1 + p_j / K_j)^n_j - 1

            each row corresponds to a reaction in the model
            each column corresponds to another set of concentrations (assuming
            lnC is a matrix)
        """
        assert lnC.shape[0] == self.Nc
        D = np.matrix(np.zeros((self.Nr, lnC.shape[1])))
        for k in range(lnC.shape[1]):
            X_k = np.log(np.exp(np.tile(lnC[:, k], (1, self.Nr))) / self.KMM + 1.0)
            ln_1_plus_S = np.matrix(np.diag(self.S_subs.T * X_k)).T
            ln_1_plus_P = np.matrix(np.diag(self.S_prod.T * X_k)).T
            D[:, k] = np.exp(ln_1_plus_S) + np.exp(ln_1_plus_P) - 1.0
        return D

    def _EtaKinetic(self, lnC):
        """
            the kinetic part of ECF3 and ECF4
        """
        return self._D_S(lnC) / self.D(lnC)

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

    def _ECF1(self, lnC):
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

    def _ECF2(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector

            Returns:
                The thermodynamic-only Enzyme Cost Function.
                Gives the predicted enzyme concentrations in [M].
        """
        assert lnC.shape == (self.Nc, 1)
        ecf2 = np.multiply(self._ECF1(lnC), 1.0/self._EtaThermodynamic(lnC))

        # fix the "fake" values that were given in ECF2 to infeasible reactions
        ecf2[ecf2 < 0] = np.nan

        return ecf2

    def _ECF3(self, lnC):
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
        return np.multiply(self._ECF2(lnC), 1.0/self._EtaKinetic(lnC))

    def _ECF4(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector

            Returns:
                The full Enzyme Cost Function, i.e. with kinetic, thermodynamic
                and allosteric data.
                Gives the predicted enzyme concentrations in [M].
        """
        assert lnC.shape == (self.Nc, 1)
        return np.multiply(self._ECF3(lnC), 1.0/self._EtaAllosteric(lnC))

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
        cap = self._ECF1(lnC)                  # capacity
        trm = 1.0/self._EtaThermodynamic(lnC)  # thermodynamics
        kin = 1.0/self._EtaKinetic(lnC)        # kinetics
        alo = 1.0/self._EtaAllosteric(lnC)     # allostery
        return np.hstack([cap, trm, kin, alo])

    def GetVolumes(self, lnC):
        """
            Arguments:
                A single metabolite ln-concentration vector

            Returns:
                Two arrays containing the enzyme volumes and
                metabolite volumes (at the provided point)
        """
        assert lnC.shape == (self.Nc, 1)
        enz_conc = self.ECF(lnC)
        met_conc = np.exp(lnC)
        enz_vol = np.multiply(enz_conc, self.mw_enz)
        met_vol = np.multiply(met_conc, self.mw_met)
        return enz_vol, met_vol

    def GetFluxes(self, lnC, E):
        assert len(lnC.shape) == 2
        assert lnC.shape[0] == self.Nc
        assert E.shape == (self.Nr, 1)

        v = np.tile(self.GetVmax(E), (1, lnC.shape[1]))
        v = np.multiply(v, self._EtaThermodynamic(lnC))
        v = np.multiply(v, self._EtaKinetic(lnC))
        v = np.multiply(v, self._EtaAllosteric(lnC))
        return v

    def ECM(self, lnC0=None, n_iter=10):
        """
            Use convex optimization to find the y with the minimal total
            enzyme cost per flux, i.e. sum(ECF(lnC))
        """

        def optfun(lnC, regularization=1e-2):
            """
                regularization function:
                    d      = x - 0.5 * (x_min + x_max)
                    lambda = median(enzyme cost weights)
                    reg    = 0.01 * lambda * 0.5 * (d.T * d)
            """
            lnC = CastToColumnVector(lnC)

            # if some reaction is not feasible, give a large penalty
            # proportional to the negative driving force.
            minimal_df = self._DrivingForce(lnC).min()
            if minimal_df <= 0:
                return 1e20 * abs(minimal_df)

            enz_conc = self.ECF(lnC)
            met_conc = np.exp(lnC)

            e = np.dot(enz_conc.T, self.mw_enz)
            m = np.dot(met_conc.T, self.mw_met)
            if np.isnan(e) or e <= 0:
                raise Exception('ECF returns NaN although all reactions are feasible')

            if self.regularization == None or self.regularization.lower() == 'none':
                return e
            elif self.regularization.lower() == 'volume':
                return e + METABOLITE_WEIGHT_CORRECTION_FACTOR * m
            elif self.regularization.lower() == 'quadratic':
                d = lnC - 0.5*(lnC.min() + lnC.max())
                return e + QUAD_REGULARIZATION_COEFF * 0.5 * float(d.T * d)
            else:
                raise Exception('Unknown regularization: ' + self.regularization)

        assert lnC0.shape == (self.Nc, 1)

        bounds = list(zip(self.lnC_bounds[:,0].flat, self.lnC_bounds[:,1].flat))

        min_res = np.inf
        lnC_min = None
        for i in range(n_iter):
            lnC0_rand = np.multiply(lnC0, 1.0 + 0.1*np.random.rand(lnC0.shape[0], 1))
            r = minimize(optfun, x0=lnC0_rand, bounds=bounds, method='SLSQP')

            if not r.success:
                continue

            res = optfun(r.x)[0, 0]
            if res < min_res:
                min_res = res
                print('%f' %  res)
                lnC_min = np.matrix(r.x).T

        return lnC_min

    def MDF(self):
        """
            Find an initial point (x0) for the optimization using MDF.
        """
        p = Pathway(self.S, self.flux, self.dG0, self.lnC_bounds)
        return p.FindMDF()

