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
from tablib.dictionary.sbtab_dict import SBtabDict
pandas.options.display.mpl_style = 'default'

exp_name = 'ecm_ecoli_aerobic'

html = HtmlWriter('res/%s.html' % exp_name)

sbtab_fpath = os.path.expanduser('~/git/enzyme-cost/data/%s.tsv' % exp_name)
sqlite_fpath = os.path.expanduser('~/git/enzyme-cost/res/%s.sqlite' % exp_name)

if not os.path.exists(sqlite_fpath):
    _sbtab_dict = SBtabDict.FromSBtab(sbtab_fpath)
    comm = sqlite3.connect(sqlite_fpath)
    _sbtab_dict.SBtab2SQL(comm)
    comm.close()

sbtab_dict = SBtabDict.FromSQLite(sqlite_fpath)
model = ECMmodel(sbtab_dict)
model.WriteMatFile('res/%s.mat' % exp_name)

lnC_MDF = model.MDF()
lnC_ECM = model.ECM()

fig1 = plt.figure(figsize=(14, 5))

ax_MDF = fig1.add_subplot(1, 2, 1)
model.PlotEnzymeCosts(lnC_MDF, ax_MDF)
ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeCosts(lnC_ECM, ax_ECM)
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