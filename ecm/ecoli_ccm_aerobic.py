# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""
import sys
import csv
from model import ECMmodel
from html_writer import HtmlWriter
import os
import sqlite3
import matplotlib.pyplot as plt
import logging
sys.path.append(os.path.expanduser('~/git/SBtab/python'))
from sqlite_interface.sbtab_dict import SBtabDict
import seaborn as sns
from errors import ThermodynamicallyInfeasibleError
import numpy as np
sns.set()
sns.axes_style("darkgrid")

l = logging.getLogger()
l.setLevel(logging.INFO)

exp_name = 'ecoli_ccm_aerobic_ProteinComposition_haverkorn'
DATA_DIR = os.path.expanduser('~/git/enzyme-cost/data')
RES_DIR = os.path.expanduser('~/git/enzyme-cost/res')
model_sbtab_fpath = os.path.join(DATA_DIR, '%s_ModelData.csv' % exp_name)
verify_sbtab_fpath = os.path.join(DATA_DIR, '%s_ValidationData.csv' % exp_name)
sqlite_fpath = os.path.join(RES_DIR, '%s.sqlite' % exp_name)
mat_fpath = os.path.join(RES_DIR, '%s.mat' % exp_name)
html_fpath = os.path.join(RES_DIR, '%s.html' % exp_name)

html = HtmlWriter(html_fpath)

logging.info('Converting input data from SBtab to SQLite')
_model_sbtab_dict = SBtabDict.FromSBtab(model_sbtab_fpath)
_verify_sbtab_dict = SBtabDict.FromSBtab(verify_sbtab_fpath)
comm = sqlite3.connect(sqlite_fpath)
_model_sbtab_dict.SBtab2SQL(comm, append=False)
_verify_sbtab_dict.SBtab2SQL(comm, append=False)
comm.commit()
comm.close()

logging.info('Reading data from SQLite database: ' + sqlite_fpath)
sbtab_dict = SBtabDict.FromSQLite(sqlite_fpath)
logging.info('Creating an ECM model using the data')
model = ECMmodel(sbtab_dict, dG0_source='keq_table')
logging.info('Exporting data to .mat file: ' + mat_fpath)
model.WriteMatFile(mat_fpath)

logging.info('Solving MDF problem')
try:
    lnC_MDF = model.MDF()
except ThermodynamicallyInfeasibleError as e:
    logging.error('The pathway is not feasible under the given constraints')
    sys.exit(-1)

USE_WOLF_CONCENTRATIONS = False
if USE_WOLF_CONCENTRATIONS:
    # load lnC data from Wolf's result matrix
    with open(os.path.join(DATA_DIR, 'metabolite_data_from_wolf.csv')) as fp:
        csv_reader = csv.DictReader(fp)
        cid2conc = {}
        for rdict in csv_reader:
            cid = rdict['!Compound:Identifiers:kegg.compound']
            cid2conc[cid] = 1e-3*float(rdict['ecf3sp'])

        cid2conc['C00001'] = 1

    lnC_ECM = np.log(np.matrix(map(cid2conc.get, model.kegg_model.cids), dtype=float).T)
else:
    # solve the ECM problem using convex optimization
    logging.info('Solving ECM problem')
    lnC_ECM = model.ECM(n_iter=15)

#fig4 = plt.figure(figsize=(14, 5))
#ax_a = fig4.add_subplot(1, 2, 1)
#model.PlotEnzymeCosts(lnC_ecf3sp, ax_a, plot_measured=True)
#ax_a.set_title(r'ecf3sp from Wolf')
#ax_b = fig4.add_subplot(1, 2, 2, sharey=ax_a)
#model.PlotEnzymeCosts(lnC_ECM, ax_b, plot_measured=True)
#ax_b.set_title(r'eCF3(1SP)')
#fig4.show()

fig1 = plt.figure(figsize=(14, 5))

ax_MDF = fig1.add_subplot(1, 2, 1)
model.PlotEnzymeCosts(lnC_MDF, ax_MDF, plot_measured=True)
ax_MDF.set_title(r'MDF')
ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeCosts(lnC_ECM, ax_ECM, plot_measured=True)
ax_ECM.set_title(r'eCF3(1SP)')
fig1.show()

fig2 = plt.figure(figsize=(14, 6))
fig2.suptitle('Valdiation')
ax_enz = fig2.add_subplot(1, 2, 1, xscale='log', yscale='log')
ax_enz.set_title('eCF3(1SP) Enzyme concentrations')
ax_met = fig2.add_subplot(1, 2, 2, xscale='log', yscale='log')
ax_met.set_title('eCF3(1SP) Metabolite concentrations')

model.ValidateMetaboliteConcentrations(lnC_ECM, ax_met)
model.ValidateEnzymeConcentrations(lnC_ECM, ax_enz)
fig2.show()

html.write('<p>\n')
html.write('<b>Experiment name:</b> %s</br>\n' % exp_name)
html.embed_matplotlib_figure(fig1)
html.write('</p><p>\n')
html.embed_matplotlib_figure(fig2)
html.write('</p><p>\n')
model.WriteHtmlTables(lnC_ECM, html)
html.write('</p>\n')

html.close()
