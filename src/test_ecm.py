# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""

from ecm import ecm_model

fpath = '/home/noore/git/enzyme-cost/data/ecm_karl.tsv'

sbtab_dict = ecm_model.SBtabDict.FromSBtab(fpath)
print sbtab_dict.reactions
print sbtab_dict.met2kegg
#cc = ComponentContribution.init()