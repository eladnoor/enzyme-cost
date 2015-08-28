# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""
import sys
from model import ECMmodel
from html_writer import HtmlWriter
import os
import sqlite3
import matplotlib.pyplot as plt
import logging
from tablib.dictionary.sbtab_dict import SBtabDict
import seaborn
seaborn.set()

l = logging.getLogger()
l.setLevel(logging.INFO)

exp_name = 'ecoli_ccm_aerobic_ProteinComposition_haverkorn'
model_sbtab_fpath = os.path.expanduser('~/git/enzyme-cost/data/%s_ModelData.csv' % exp_name)
verify_sbtab_fpath = os.path.expanduser('~/git/enzyme-cost/data/%s_ValidationData.csv' % exp_name)
sqlite_fpath = os.path.expanduser('~/git/enzyme-cost/res/%s.sqlite' % exp_name)
mat_fpath = os.path.expanduser('~/git/enzyme-cost/res/%s.mat' % exp_name)
html_fpath = os.path.expanduser('~/git/enzyme-cost/res/%s.html' % exp_name)

html = HtmlWriter(html_fpath)

#if not os.path.exists(sqlite_fpath):
if True: # always override the SQL database
    logging.info('Converting input data from SBtab to SQLite')
    _model_sbtab_dict = SBtabDict.FromSBtab(model_sbtab_fpath)
    _verify_sbtab_dict = SBtabDict.FromSBtab(verify_sbtab_fpath)

    os.remove(sqlite_fpath)
    comm = sqlite3.connect(sqlite_fpath)
    _model_sbtab_dict.SBtab2SQL(comm)
    _verify_sbtab_dict.SBtab2SQL(comm)

    comm.close()

logging.info('Reading data from SQLite database: ' + sqlite_fpath)
sbtab_dict = SBtabDict.FromSQLite(sqlite_fpath)
logging.info('Creating an ECM model using the data')
model = ECMmodel(sbtab_dict, calculate_dG0_using_CC=True)
logging.info('Exporting data to .mat file: ' + mat_fpath)
model.WriteMatFile(mat_fpath)

logging.info('Solving MDF problem')
lnC_MDF = model.MDF()
logging.info('Solving ECM problem')
lnC_ECM = model.ECM()

fig1 = plt.figure(figsize=(14, 5))

ax_MDF = fig1.add_subplot(1, 2, 1)
model.PlotEnzymeCosts(lnC_MDF, ax_MDF, plot_measured=True)
ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeCosts(lnC_ECM, ax_ECM, plot_measured=True)
fig1.show()

fig2 = plt.figure(figsize=(6, 6))
fig2.suptitle('Metabolite Concentrations')
ax = fig2.add_subplot(1, 1, 1, xscale='log', yscale='log')
model.ValidateMetaboliteConcentrations(lnC_ECM, ax)
fig2.show()

fig3 = plt.figure(figsize=(6, 6))
fig3.suptitle('Enzyme Concentrations')
ax = fig3.add_subplot(1, 1, 1, xscale='log', yscale='log')
model.ValidateEnzymeConcentrations(lnC_ECM, ax)
fig3.show()

html.write('<p>\n')
html.write('<b>Experiment name:</b> %s</br>\n' % exp_name)
html.embed_matplotlib_figure(fig1)
html.write('</p><p>\n')
html.embed_matplotlib_figure(fig2)
html.embed_matplotlib_figure(fig3)
html.write('</p><p>\n')
model.WriteHtmlTables(lnC_ECM, html)
html.write('</p>\n')

html.close()