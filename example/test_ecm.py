# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""

from ecm.model import ECMmodel
import os
import numpy as np
import matplotlib.pyplot as plt

fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_karl.tsv')

model = ECMmodel(fpath)

y1 = model.ecf.MDF()
y2 = model.ecf.ECM()

plt.plot(y1, y2, '.')
#print '\n'.join(map(lambda (cid,x): '%s = %.2e' % (cid, np.exp(x)), zip(model.kegg_model.cids, y1.flat)))
#print '\n'.join(map(lambda (cid,x): '%s = %.2e' % (cid, np.exp(x)), zip(model.kegg_model.cids, y2.flat)))