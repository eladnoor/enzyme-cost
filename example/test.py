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

fpath = os.path.expanduser('~/git/enzyme-cost/data/ecm_karl.tsv')

model = ECMmodel(fpath)
model.WriteMatFile('res/karl.mat')

fig = plt.figure(figsize=(15, 7))

lnC_MDF = model.MDF()
ax_MDF = fig.add_subplot(1, 2, 1)
model.PlotEnzymeCosts(lnC_MDF, ax_MDF)

lnC_ECM = model.ECM()
ax_ECM = fig.add_subplot(1, 2, 2, sharey=ax_MDF)
model.PlotEnzymeCosts(lnC_ECM, ax_ECM)

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(1, 1, 1)
model.ValidateMetaboliteConcentrations(lnC_ECM, ax)

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(1, 1, 1)
model.ValidateEnzymeConcentrations(lnC_ECM, ax)