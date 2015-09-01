# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:30:45 2015

@author: eladn
"""
import numpy as np
from collections import Iterable
import re
from scipy import stats

def CastToColumnVector(v):
    """
        casts any numeric list of floats to a 2D-matrix with 1 column,
        and rows corresponding to the length of the list
        
        if the input is a NumPy array or matrix, it will be reshaped to
    """
    if type(v) in [np.ndarray, np.matrix]:
        return np.matrix(np.reshape(v, (np.prod(v.shape), 1)), dtype=float)
    if isinstance(v, Iterable):
        return np.matrix(list(v), dtype=float).T
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

def PlotCorrelation(ax, x, y, labels, mask=None, scale='log', grid=True):
    """
        scale - if 'log' indicates that the regression should be done on the
                logscale data.
    """
    x = CastToColumnVector(x)
    y = CastToColumnVector(y)
    
    if mask is None:
        mask = (x > 0) & (y > 0)
    
    ax.grid(grid)
    if scale == 'log':
        ax.set_xscale('log')
        ax.set_yscale('log')
        slope, intercept, r_value, p_value, std_err = \
           stats.linregress(np.log(x[mask]).flat, np.log(y[mask]).flat)
    else:
        ax.set_xscale('linear')
        ax.set_yscale('linear')
        slope, intercept, r_value, p_value, std_err = \
            stats.linregress(x[mask].flat, y[mask].flat)
    ax.plot(x[mask], y[mask], '.', markersize=15, color='red', alpha=0.5)
    ax.plot(x[~mask], y[~mask], '.', markersize=15, color='blue', alpha=0.5)
    
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_xlim(min(xmin, ymin), max(xmax, ymax))
    ax.set_ylim(min(xmin, ymin), max(xmax, ymax))
    ax.plot([0, 1], [0, 1], ':', color='black', alpha=0.4, transform=ax.transAxes)
        
    ax.text(0.05, 0.9, '$r$ = %.2f\n(p = %.1e)' % (r_value, p_value),
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes,
            color='black', fontsize=10,
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':10}) 
    
    for l, x_i, y_i, m in zip(labels, x, y, mask):
        if m:
            ax.text(x_i, y_i, l, alpha=1.0)
        elif np.isfinite(x_i) and np.isfinite(y_i):
            if scale == 'linear' or (x_i > 0 and y_i > 0):
                ax.text(x_i, y_i, l, alpha=0.4)

if __name__ == '__main__':
    x = np.eye(3)
    print CastToColumnVector(x)
    print CastToColumnVector([4,3,1,2])