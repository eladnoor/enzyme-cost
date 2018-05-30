#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 17:22:11 2017

@author: noore
"""

import unittest
import warnings
import numpy as np
import pandas as pd
from ecm.model import ECMmodel
from ecm.cost_function import EnzymeCostFunction
from ecm.optimized_bottleneck_driving_force import Pathway

def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            test_func(self, *args, **kwargs)
    return do_test

class TestReactionParsing(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestReactionParsing, self).__init__(*args, **kwargs)
        R = 8.31e-3
        DEFAULT_TEMP = 298.15  # K
        RT = R * DEFAULT_TEMP
        
        Nr = 3
        Nc = 4
        
        self.S = np.matrix(np.zeros((Nc, Nr)))
        
        self.S[0, 0] = -1
        self.S[1, 0] = 1
        self.S[2, 0] = 1
        self.S[1, 1] = -1
        self.S[2, 1] = 1
        self.S[2, 2] = -1
        self.S[3, 2] = 1
        
        self.v          = np.matrix([1.0, 1.0, 2.0]).T
        self.kcat       = np.matrix([1.0, 1.0, 1.0]).T
        self.dGm_r      = np.matrix([-3.0, -2.0, -3.0]).T
        self.dG0_r      = self.dGm_r - RT * self.S.T * np.matrix(np.ones((Nc, 1))) * np.log(1e-3)
        self.K_M        = np.matrix(np.ones(self.S.shape))
        self.K_M[self.S < 0] = 9e-2
        self.K_M[self.S > 0] = 1e-2
        
        self.lnC_bounds = np.log(np.matrix([[1e-9]*Nc, [1e-1]*Nc]).T)
        
        self.A_act      = np.matrix(np.zeros(self.S.shape))
        self.A_inh      = np.matrix(np.zeros(self.S.shape))
        self.K_act      = np.matrix(np.ones(self.S.shape))
        self.K_inh      = np.matrix(np.ones(self.S.shape))

    @ignore_warnings
    def test_toy_mdf(self):
        
        toy_ecf = Pathway(self.S, self.v, self.dG0_r, self.lnC_bounds)
        mdf, params = toy_ecf.FindMDF()
        self.assertAlmostEqual(mdf, 9.169, 3)
        
    @ignore_warnings
    def test_toy_ecm(self):
        
        toy_ecf = EnzymeCostFunction(self.S, self.v, self.kcat, self.dG0_r,
                                     self.K_M, self.lnC_bounds, None, None,
                                     self.A_act, self.A_inh,
                                     self.K_act, self.K_inh)
        lnC = toy_ecf.ECM(n_iter=20)
        concs = np.exp(lnC)
        costs = toy_ecf.ECF(lnC)
        self.assertAlmostEqual(concs[0,0], 0.1, 2)
        self.assertAlmostEqual(concs[1,0], 0.015, 3)
        self.assertAlmostEqual(concs[2,0], 0.0067, 4)
        self.assertAlmostEqual(concs[3,0], 1e-9, 8)

        self.assertAlmostEqual(costs[0,0], 3.93, 2)
        self.assertAlmostEqual(costs[1,0], 3.08, 2)
        self.assertAlmostEqual(costs[2,0], 5.22, 2)

    @ignore_warnings
    def test_ecm(self):
        df_names = ECMmodel.DATAFRAME_NAMES
        df_dict = {n : pd.read_csv('tests/test_%s.csv' % n) for n in df_names}
        ecm_model = ECMmodel(df_dict, flux_unit='mM/s', bound_unit='M')
        lnC = ecm_model.ECM()
        score = ecm_model.ecf.ECF(lnC).sum()
        self.assertAlmostEqual(score, 0.01514, 5)
        
if __name__ == '__main__':
    unittest.main()
