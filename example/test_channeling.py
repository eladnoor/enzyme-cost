# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""
from ecm.model import ECMmodel
from ecm.html_writer import HtmlWriter
import os
import sqlite3
import matplotlib.pyplot as plt
import pandas
import logging
import sys
sys.path.append(os.path.expanduser('~/git/SBtab/python'))
from sqlite_interface.sbtab_dict import SBtabDict
pandas.options.display.mpl_style = 'default'

l = logging.getLogger()
l.setLevel(logging.INFO)

exp_name = 'ecm_ecoli_aerobic_channeling'
sbtab_fpath = os.path.expanduser('~/git/enzyme-cost/data/%s.tsv' % exp_name)
sqlite_fpath = os.path.expanduser('~/git/enzyme-cost/res/%s.sqlite' % exp_name)
mat_fpath = os.path.expanduser('~/git/enzyme-cost/res/%s.mat' % exp_name)
html_fpath = os.path.expanduser('~/git/enzyme-cost/res/%s.html' % exp_name)

html = HtmlWriter(html_fpath)

if not os.path.exists(sqlite_fpath):
    logging.info('Converting input data from SBtab to SQLite')
    _sbtab_dict = SBtabDict.FromSBtab(sbtab_fpath)
    comm = sqlite3.connect(sqlite_fpath)
    _sbtab_dict.SBtab2SQL(comm)
    comm.close()

logging.info('Reading data from SQLite database: ' + sqlite_fpath)
sbtab_dict = SBtabDict.FromSQLite(sqlite_fpath)
logging.info('Creating an ECM model using the data')
model = ECMmodel(sbtab_dict, dG0_source='component_contribution')
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