# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:38:13 2015

@author: noore
"""

from SBtabTools import oneOrMany
from SBtab import SBtabTable
from tablibIO import loadTSV

class ECM_Model(object):
    
    def __init__(self, sbtab):
        self.sbtab = sbtab
        pass
    
    @staticmethod
    def FromSBtab(fpath):
        spreadsheet_file = loadTSV(fname, False)
        m = oneOrMany(spreadsheet_file)
        sbtab_list = [SBtabTable(dset, fname) for dset in m]
        return sbtab_list
        
if __name__ == '__main__':
    fname = '/home/noore/git/enzyme-cost/data/ecm_karl.tsv'
    m = ECM_Model.FromSBtab(fname)
    