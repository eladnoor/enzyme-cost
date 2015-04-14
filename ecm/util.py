# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:30:45 2015

@author: eladn
"""
import numpy as np
import types
from collections import Iterable

def CastToColumnVector(v):
    """
        casts any numeric list of floats to a 2D-matrix with 1 column,
        and rows corresponding to the length of the list
        
        if the input is a NumPy array or matrix, it will be reshaped to
    """
    if type(v) in [np.ndarray, np.matrix]:
        return np.matrix(np.reshape(v, (np.prod(v.shape), 1)))
    if isinstance(v, Iterable):
        return np.matrix(iter(v)).T
    else:
        raise ValueError('Can only cast lists or numpy arrays, not ' + str(type(v)))


if __name__ == '__main__':
    x = np.eye(3)
    print CastToColumnVector(x)
    print CastToColumnVector([4,3,1,2])