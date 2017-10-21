# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:30:45 2015

@author: eladn
"""
import numpy as np
from collections import Iterable
from scipy import stats
import urllib

RT = 8.31e-3 * 298.15
CELL_VOL_PER_DW = 2.7e-3 # L/gCDW [Winkler and Wilson, 1966, http://www.jbc.org/content/241/10/2200.full.pdf+html]
ECF_DEFAULTS = {
    'version'       : 3,
    'dG0_source'    : 'keq_table', # options are: 'keq_table', 'dG0r_table', or 'component_contribution'
    'kcat_source'   : 'gmean',     # options are: 'fwd' or 'gmean'
    'denominator'   : 'CM',        # options are: 'S', 'SP', '1S', '1SP', or 'CM'
    'regularization': 'volume',    # options are: None, 'volume', or 'quadratic'
    }
str2bool = lambda x : x not in [0, 'false', 'False']

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
        logx = np.log10(x[mask])
        logy = np.log10(y[mask])
        slope, intercept, r_value, p_value, std_err = \
           stats.linregress(logx.flat, logy.flat)
        rmse = np.sqrt( np.power(logx - logy, 2).mean() )
    else:
        ax.set_xscale('linear')
        ax.set_yscale('linear')
        slope, intercept, r_value, p_value, std_err = \
            stats.linregress(x[mask].flat, y[mask].flat)
        rmse = np.sqrt( np.power(x[mask] - y[mask], 2).mean() )
    ax.plot(x[mask], y[mask], '.', markersize=15, color='red', alpha=0.5)
    ax.plot(x[~mask], y[~mask], '.', markersize=15, color='blue', alpha=0.5)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_xlim(min(xmin, ymin), max(xmax, ymax))
    ax.set_ylim(min(xmin, ymin), max(xmax, ymax))
    ax.plot([0, 1], [0, 1], ':', color='black', alpha=0.4, transform=ax.transAxes)

    boxtext = 'RMSE = %.2f\n$r^2$ = %.2f\n(p = %.1e)' % \
              (rmse, r_value**2, p_value)

    ax.text(0.05, 0.95, boxtext,
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes,
            color='black', fontsize=10,
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    for l, x_i, y_i, m in zip(labels, x, y, mask):
        if m:
            ax.text(x_i, y_i, l, alpha=1.0)
        elif np.isfinite(x_i) and np.isfinite(y_i):
            if scale == 'linear' or (x_i > 0 and y_i > 0):
                ax.text(x_i, y_i, l, alpha=0.4)

def CompoundID2MolecularWeight(compound_id):
    import openbabel
    s_mol = urllib.urlopen('http://rest.kegg.jp/get/cpd:%s/mol' % compound_id).read()
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('mol', 'inchi')
    conv.AddOption("F", conv.OUTOPTIONS)
    conv.AddOption("T", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s_mol)):
        return None
    else:
        return obmol.GetMolWt()

if __name__ == '__main__':
    pass