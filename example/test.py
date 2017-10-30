# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""
from ecm.model import ECMmodel
from util.SBtab.SBtabTools import SBtabDict
import os
import matplotlib.pyplot as plt
import logging
import inspect

l = logging.getLogger()
l.setLevel(logging.INFO)

SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
DATA_DIR = os.path.join(os.path.split(SCRIPT_DIR)[0], 'data')
RESULT_DIR = os.path.join(os.path.split(SCRIPT_DIR)[0], 'res')

#exp_name = 'ecoli_ccm_aerobic_ProteinComposition_haverkorn'
exp_name = 'ecoli_ccm_aerobic_channeling'
modeldata_fname = os.path.join(DATA_DIR, '%s_ModelData.tsv' % exp_name)
validationdata_fname = os.path.join(DATA_DIR, '%s_ValidationData.tsv' % exp_name)

logging.info('Reading SBtab files')
modeldata_sbtabs = SBtabDict.FromSBtab(modeldata_fname)
validationdata_sbtabs = SBtabDict.FromSBtab(validationdata_fname)

logging.info('Creating an ECM model using the data')
model = ECMmodel(modeldata_sbtabs, validationdata_sbtabs)

logging.info('Solving MDF problem')
lnC_MDF = model.MDF()
logging.info('Solving ECM problem')
lnC_ECM = model.ECM(n_iter=5)

fig1 = plt.figure(figsize=(14, 5))
ax_MDF = fig1.add_subplot(1, 2, 1)
model.PlotEnzymeDemandBreakdown(lnC_MDF, ax_MDF, plot_measured=True)
ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeDemandBreakdown(lnC_ECM, ax_ECM, plot_measured=True)

fig2 = plt.figure(figsize=(6, 6))
fig2.suptitle('Metabolite Concentrations')
ax = fig2.add_subplot(1, 1, 1, xscale='log', yscale='log')
model.ValidateMetaboliteConcentrations(lnC_ECM, ax)

fig3 = plt.figure(figsize=(6, 6))
fig3.suptitle('Enzyme Concentrations')
ax = fig3.add_subplot(1, 1, 1, xscale='log', yscale='log')
model.ValidateEnzymeConcentrations(lnC_ECM, ax)

#%%
fig4 = plt.figure(figsize=(6, 3))
ax_ECM = fig4.add_subplot(1, 1, 1, sharey=ax_MDF)
model.PlotEnzymeDemandBreakdown(lnC_ECM, ax_ECM, plot_measured=True)
fig4.savefig(os.path.join(RESULT_DIR, 'breakdown.svg'))

fig5 = plt.figure(figsize=(5, 5))
ax = fig5.add_subplot(1, 1, 1)
#model.PlotVolumesPie(lnC_ECM, ax)
vols, labels, colors = model._GetVolumeDataForPlotting(lnC_ECM)
ax.pie(vols, labels=labels, colors=colors)
ax.set_title('total weight = %.2g [g/L]' % sum(vols))

fig5.savefig(os.path.join(RESULT_DIR, 'pie.svg'))
