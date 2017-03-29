# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""
import sys
from model import ECMmodel
import os
import matplotlib.pyplot as plt
import logging
sys.path.append(os.path.expanduser('~/git/SBtab/python'))
from sqlite_interface.sbtab_dict import SBtabDict
import seaborn as sns
from errors import ThermodynamicallyInfeasibleError
sns.set()
sns.axes_style("darkgrid")

l = logging.getLogger()
l.setLevel(logging.INFO)

exp_name = 'ecoli_ccm_aerobic_ProteinComposition_haverkorn'
DATA_DIR = os.path.expanduser('~/git/enzyme-cost/data')
RES_DIR = os.path.expanduser('~/git/enzyme-cost/res')
model_sbtab_fpath = os.path.join(DATA_DIR, '%s_ModelData.csv' % exp_name)
validate_sbtab_fpath = os.path.join(DATA_DIR, '%s_ValidationData.csv' % exp_name)
mat_fpath = os.path.join(RES_DIR, '%s.mat' % exp_name)

logging.info('Converting input data from SBtab to SQLite')
_model_sbtab_dict = SBtabDict.FromSBtab(model_sbtab_fpath)
_validate_sbtab_dict = SBtabDict.FromSBtab(validate_sbtab_fpath)

logging.info('Creating an ECM model using the data')
model = ECMmodel(_model_sbtab_dict, _validate_sbtab_dict,
                 ecf_version=3,
                 denom_version='CM',
                 regularization='volume',
                 dG0_source='keq_table',
                 kcat_source='gmean')
ecf_title = 'ECF3(CM)'

logging.info('Exporting data to .mat file: ' + mat_fpath)
model.WriteMatFile(mat_fpath)

logging.info('Solving MDF problem')
try:
    lnC_MDF = model.MDF()
except ThermodynamicallyInfeasibleError as e:
    logging.error('The pathway is not feasible under the given constraints')
    sys.exit(-1)

# solve the ECM problem using convex optimization
logging.info('Solving ECM problem')
lnC_ECM = model.ECM(n_iter=15)

fig1 = plt.figure(figsize=(12, 4))

ax_MDF = fig1.add_subplot(1, 2, 1)
model.PlotEnzymeCosts(lnC_MDF, ax_MDF, plot_measured=True)
ax_MDF.set_title(r'MDF')
ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeCosts(lnC_ECM, ax_ECM, plot_measured=True)
ax_ECM.set_title(ecf_title)
fig1.savefig(os.path.join(RES_DIR, '%s_bar.pdf' % exp_name))

fig2 = plt.figure(figsize=(12, 5))
fig2.suptitle('Valdiation')
ax_enz = fig2.add_subplot(1, 2, 1, xscale='log', yscale='log')
ax_enz.set_title('%s Enzyme concentrations' % ecf_title)
ax_met = fig2.add_subplot(1, 2, 2, xscale='log', yscale='log')
ax_met.set_title('%s Metabolite concentrations' % ecf_title)

model.ValidateMetaboliteConcentrations(lnC_ECM, ax_met)
model.ValidateEnzymeConcentrations(lnC_ECM, ax_enz)
fig2.savefig(os.path.join(RES_DIR, '%s_corr.pdf' % exp_name))
