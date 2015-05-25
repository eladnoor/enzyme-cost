# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:31:55 2015

@author: noore
"""

from SBtabTools import oneOrMany
from SBtab import SBtabTable, SBtabError
from tablibIO import loadTSV
import os

fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_ecoli_aerobic.tsv')

spreadsheet_file = loadTSV(fpath, False)
m = oneOrMany(spreadsheet_file)
reaction_sbtab = SBtabTable(m[0], fpath)
