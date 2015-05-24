# -*- coding: utf-8 -*-
"""
Created on Thu May 14 17:00:18 2015

@author: eladn
"""
import os
from ecm.sbtab_dict import SBtabDict
import sqlite3

sbtab_fname = os.path.expanduser('~/git/enzyme-cost/data/ecm_ecoli_aerobic.tsv')
sqlite_fname = os.path.expanduser('~/git/enzyme-cost/res/ecm_ecoli_aerobic.sqlite')

if False:
    _sbtab_dict = SBtabDict.FromSBtab(sbtab_fname)
    comm = sqlite3.connect(sqlite_fname)
    _sbtab_dict.SBtab2SQL(comm)
    comm.close()
else:
    _sbtab_dict = SBtabDict.FromSQLite(sqlite_fname)


# TODO:
# sqlite -> SBtab export function
# API:
#     - GetKineticData(organism_name)
#     - GetExperimentalData(data_type, organism_name, experiment_id)
#     - GetBooleanGeneFormula(reaction_name)
# potentially, join tables with their accompanying table with titles
#
# convert Wolf's table to sqlite and upload to GitHub
#
# gene(orf) => complex <=> EC number => reaction ID = reaction formula