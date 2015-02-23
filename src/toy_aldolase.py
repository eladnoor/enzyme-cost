# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 11:43:14 2015

@author: noore

Description:
    create a toy diamon network model (i.e. two parallel pathways from A to D, via
    intermediates B or C)

"""

import numpy as np
from ecf import ECF, RT

import os
if not os.path.exists('../res'):
    os.mkdir('../res')


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

v      = np.matrix([1.0, 1.0, 2.0]).T
kcat   = np.matrix([1.0, 1.0, 1.0]).T
dGm_r  = np.matrix([-3.0, -2.0, -3.0]).T
dG0_r  = dGm_r - RT * S.T * np.matrix(np.ones((Nc, 1))) * np.log(1e-3)
K_subs = 9e-2*np.matrix(np.matrix(np.ones((Nc, Nr))))
K_prod = 1e-2*np.matrix(np.matrix(np.ones((Nc, Nr))))

toy_ecf = ECF(S, v, kcat, dG0_r, K_subs, K_prod)
toy_ecf.generate_pdf_report('../res/aldolase.pdf')
