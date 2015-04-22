# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:59 2015

@author: noore
"""

from ecm.model import ECMmodel
from ecm.html_writer import HtmlWriter
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas
pandas.options.display.mpl_style = 'default'

html = HtmlWriter('res/karl.html')

fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_karl.tsv')
#fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_ecoli_aerobic.tsv')

model = ECMmodel(fpath)
#model.WriteMatFile('res/karl.mat')

lnC_MDF = model.MDF()
lnC_ECM = model.ECM()

fig1 = plt.figure(figsize=(15, 7))

ax_MDF = fig1.add_subplot(1, 2, 1)
model.PlotEnzymeCosts(lnC_MDF, ax_MDF)
ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeCosts(lnC_ECM, ax_ECM)
html.embed_matplotlib_figure(fig1)
fig1.show()

fig2 = plt.figure(figsize=(10, 10))
fig2.suptitle('Metabolite Concentrations')
ax = fig2.add_subplot(1, 1, 1, xscale='log', yscale='log')
model.ValidateMetaboliteConcentrations(lnC_ECM, ax)
fig2.show()
html.embed_matplotlib_figure(fig2)

fig3 = plt.figure(figsize=(10, 10))
fig3.suptitle('Enzyme Concentrations')
ax = fig3.add_subplot(1, 1, 1, xscale='log', yscale='log')
model.ValidateEnzymeConcentrations(lnC_ECM, ax)
fig3.show()
html.embed_matplotlib_figure(fig3)

met_conc = dict(zip(model.kegg_model.cids, np.exp(lnC_ECM).flat))
enz_conc = dict(zip(model.kegg_model.rids, model.ecf.ECF(lnC_ECM).flat))

html.write_table([{'reaction':          r,
                   'flux [mM/s]':       model.rid2flux[r]*1e3,
                   'enzyme conc. [mM]': enz_conc[r]*1e3}
                    for r in model.kegg_model.rids],
                 headers=['reaction', 'flux [mM/s]', 'enzyme conc. [mM]'],
                 decimal=4)
                 
html.write_table([{'compound':          cid,
                   'concentration [mM]': met_conc[cid]*1e3}
                  for cid in model.kegg_model.cids],
                 headers=['compound', 'concentration [mM]'],
                 decimal=4)

html.close()