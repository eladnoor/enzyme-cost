# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""

from ecm.ecm_model import ECMmodel
import os

fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_karl.tsv')

ecm_model = ECMmodel(fpath)

for i, rid in enumerate(ecm_model.kegg_model.rids):
    #print rid, ecm_model.kegg_model.write_reaction_by_index(i)
    print '%50s: %.1f, %.1f, %.1f' % (ecm_model.kegg_model.write_reaction_by_index(i), 
                            ecm_model.rid2dG0_cc[rid],
                            ecm_model.rid2dG0_gibbs_energy_table[rid],
                            ecm_model.rid2dG0_rate_constant_table[rid],
                            )
    