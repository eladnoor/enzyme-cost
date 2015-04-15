# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:30:45 2015

@author: eladn
"""
import numpy as np
from collections import Iterable
import re

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

def ParseReactionFormulaSide(s):
    """ 
        Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        Ignores stoichiometry.
        
        Returns:
            The set of CIDs.
    """
    if s.strip() == "null":
        return {}
    
    compound_bag = {}
    for member in re.split('\s+\+\s+', s):
        tokens = member.split(None, 1)
        if len(tokens) == 0:
            continue
        if len(tokens) == 1:
            amount = 1
            key = member
        else:
            try:
                amount = float(tokens[0])
            except ValueError:
                raise Exception(
                    "Non-specific reaction: %s" % s)
            key = tokens[1]
            
        try:
            compound_bag[key] = compound_bag.get(key, 0) + amount
        except ValueError:
            raise Exception(
                "Non-specific reaction: %s" % s)
    
    return compound_bag

def ParseReaction(formula, arrow='<=>'):
    """ 
        Parses a two-sided formula such as: 2 FBP => DHAP + GAP
        
        Return:
            The set of substrates, products and the direction of the reaction
    """
    tokens = formula.split(arrow)
    if len(tokens) < 2:
        raise Exception('Reaction does not contain the arrow sign (%s): %s'
                        % (arrow, formula))
    if len(tokens) > 2:
        raise Exception('Reaction contains more than one arrow sign (%s): %s'
                        % (arrow, formula))
    
    left = tokens[0].strip()
    right = tokens[1].strip()
    
    sparse_reaction = {}
    for cid, count in ParseReactionFormulaSide(left).iteritems():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 

    for cid, count in ParseReactionFormulaSide(right).iteritems():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 

    return sparse_reaction

if __name__ == '__main__':
    x = np.eye(3)
    print CastToColumnVector(x)
    print CastToColumnVector([4,3,1,2])