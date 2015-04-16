# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""

from ecm.model import ECMmodel
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas
pandas.options.display.mpl_style = 'default'

#fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_karl.tsv')
fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_ecoli_aerobic.tsv')

model = ECMmodel(fpath, thermo_mode='CC')
model.WriteMatFile('res/karl.mat')

lnC_MDF = model.MDF()
lnC_ECM = model.ECM()

fig1 = plt.figure(figsize=(15, 7))

ax_MDF = fig1.add_subplot(1, 2, 1)
model.PlotEnzymeCosts(lnC_MDF, ax_MDF)
ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeCosts(lnC_ECM, ax_ECM)
fig1.show()

fig2 = plt.figure(figsize=(10, 10))
fig2.suptitle('Metabolite Concentrations')
ax = fig2.add_subplot(1, 1, 1)
model.ValidateMetaboliteConcentrations(lnC_ECM, ax)
fig2.show()

fig3 = plt.figure(figsize=(10, 10))
fig3.suptitle('Enzyme Concentrations')
ax = fig3.add_subplot(1, 1, 1)
model.ValidateEnzymeConcentrations(lnC_ECM, ax)
fig3.show()


met_conc = np.exp(lnC_ECM)

enz_conc = dict(zip(model.kegg_model.rids, model.ecf.ECF(lnC_ECM).flat))
