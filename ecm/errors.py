# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:30:38 2015

@author: noore
"""
import numpy as np

class NonSteadyStateSolutionError(Exception):
    pass

class ThermodynamicallyInfeasibleError(Exception):

    def __init__(self, lnC=None):
        if lnC is None:
            Exception.__init__(self, 'this reaction system has no feasible solutions')
        else:
            Exception.__init__(self, 'C = %s : is thermodynamically infeasible' 
                                     % str(np.exp(lnC)))
    pass
