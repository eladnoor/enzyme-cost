# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""

from ecm.model import ECMmodel
import os
import numpy as np

fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_karl.tsv')

model = ECMmodel(fpath)

y = model.ecf.MDF()
print '\n'.join(map(lambda (cid,x): '%s = %.2e' % (cid, np.exp(x)), zip(model.kegg_model.cids, y.flat)))