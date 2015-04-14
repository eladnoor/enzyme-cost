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

for func in [model.MDF, model.ECM]:
    lnC = func()
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(311, yscale='log')
    ax.plot(np.exp(lnC), 'o')
    ax.set_xlabel('metabolite')
    ax.set_ylabel('concentration [M]')
    ax = fig.add_subplot(312)
    ax.plot(model.ecf._DrivingForce(lnC), 'o')
    ax.set_xlabel('enzyme')
    ax.set_ylabel('driving force [kJ/mol]')
    ax = fig.add_subplot(313, yscale='linear')
    model.PlotEnzymeCosts(lnC, ax)
