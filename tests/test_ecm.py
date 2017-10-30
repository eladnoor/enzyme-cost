#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 17:22:11 2017

@author: noore
"""

import unittest
import numpy as np
from ecm.cost_function import EnzymeCostFunction

class TestReactionParsing(unittest.TestCase):
    
    def test_rump_pulp(self):
        R = 8.31e-3
        DEFAULT_TEMP = 298.15  # K
        RT = R * DEFAULT_TEMP
        
        Nr = 3
        Nc = 4
        
        S = np.matrix(np.zeros((Nc, Nr)))
        
        S[0, 0] = -1
        S[1, 0] = 1
        S[2, 0] = 1
        S[1, 1] = -1
        S[2, 1] = 1
        S[2, 2] = -1
        S[3, 2] = 1
        
        v          = np.matrix([1.0, 1.0, 2.0]).T
        kcat       = np.matrix([1.0, 1.0, 1.0]).T
        dGm_r      = np.matrix([-3.0, -2.0, -3.0]).T
        dG0_r      = dGm_r - RT * S.T * np.matrix(np.ones((Nc, 1))) * np.log(1e-3)
        K_M        = np.matrix(np.ones(S.shape))
        K_M[S < 0] = 9e-2
        K_M[S > 0] = 1e-2
        
        lnC_bounds = np.log(np.matrix([[1e-9]*Nc, [1e-1]*Nc]).T)
        
        A_act      = np.matrix(np.zeros(S.shape))
        A_inh      = np.matrix(np.zeros(S.shape))
        K_act      = np.matrix(np.ones(S.shape))
        K_inh      = np.matrix(np.ones(S.shape))
        
        toy_ecf = EnzymeCostFunction(S, v, kcat, dG0_r, K_M, lnC_bounds, None, None,
                                     A_act, A_inh, K_act, K_inh)
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
        
if __name__ == '__main__':
    unittest.main()
