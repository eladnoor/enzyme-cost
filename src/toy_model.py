# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 12:15:06 2015

@author: noore

Description:
    Generate a toy reaction network model with N intermediate metabolites
    and N+1 reactions.
    
"""

from ecf import ECF
import numpy as np

import os
if not os.path.exists('../res'):
    os.mkdir('../res')

N = 2

S = np.matrix(np.zeros((N + 2, N + 1)))
for i in xrange(N + 1):
    S[i, i] = -1.0
    S[i+1, i] = 1.0
v        = np.matrix(np.ones((N + 1, 1)))
kcat     = np.matrix(np.ones((N + 1, 1)))
dG0      = -5.0 * np.matrix(np.ones((N + 1, 1)))
K_M      = np.matrix(np.ones(S.shape))
K_M[S < 0] = 9e-2
K_M[S > 0] = 1e-2

ecf = ECF(S, v, kcat, dG0, K_M)

ecf.generate_pdf_report('../res/toy.pdf')
