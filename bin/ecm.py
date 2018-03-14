#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 17:44:09 2017

@author: noore
"""
import argparse
import os
import logging
import inspect
import sys
sys.path.append(os.path.expanduser('~/git/SBtab/python'))
sys.path.append(os.path.expanduser('~/git/enzyme-cost'))
from sqlite_interface.sbtab_dict import SBtabDict
from ecm.model import ECMmodel
import json
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
BASE_DIR = os.path.split(SCRIPT_DIR)[0]

def MakeParser():
    parser = argparse.ArgumentParser(
        description=('Run Enzyme Cost Minimization (ECM)'))
    parser.add_argument('config_fname', type=str,
                        help='Configuration filename')
    return parser

###############################################################################
parser = MakeParser()
args = parser.parse_args()
with open(args.config_fname, 'r') as fp:
    config = json.load(fp)

input_sbtab_fname = \
    os.path.expanduser(config['IO']['model_sbtab_fname'])
validation_sbtab_fname = \
    os.path.expanduser(config['IO'].get('validation_sbtab_fname', ''))
output_prefix = \
    os.path.expanduser(config['IO']['output_prefix'])

ecf_params = config['ECF']

logging.getLogger().setLevel(logging.INFO)

logging.info('Reading SBtab file: ' + input_sbtab_fname)
modeldata_sbtabs = SBtabDict.FromSBtab(input_sbtab_fname)

logging.info('Creating an ECM model using the data')

if validation_sbtab_fname:
    validationdata_sbtabs = SBtabDict.FromSBtab(validation_sbtab_fname)
else:
    validationdata_sbtabs = None
model = ECMmodel(modeldata_sbtabs,
                 validate_sbtab=validationdata_sbtabs,
                 ecf_params=ecf_params)

logging.info('Solving MDF problem')
lnC_MDF = model.MDF()
logging.info('Solving ECM problem')
lnC_ECM = model.ECM(n_iter=5)

if config['IO']['generate_result_sbtab']:
    res_sbtab = model.ToSBtab(lnC_ECM,
                              output_prefix,
                              document_name='E. coli central carbon metabolism - ECM result')

if config['IO']['generate_result_figures']:
    fig1 = plt.figure(figsize=(14, 5))
    ax_MDF = fig1.add_subplot(1, 2, 1)
    model.PlotEnzymeDemandBreakdown(lnC_MDF, ax_MDF, plot_measured=True)
    ax_MDF.set_title('MDF results')
    ax_ECM = fig1.add_subplot(1, 2, 2, sharey=ax_MDF)
    model.PlotEnzymeDemandBreakdown(lnC_ECM, ax_ECM, plot_measured=True)
    ax_ECM.set_title('ECF results')
    fig1.savefig(output_prefix + 'enzyme_demand.svg')
    
    if validation_sbtab_fname:
        fig2 = plt.figure(figsize=(6, 6))
        fig2.suptitle('Metabolite Concentrations')
        ax = fig2.add_subplot(1, 1, 1, xscale='log', yscale='log')
        model.ValidateMetaboliteConcentrations(lnC_ECM, ax)
        fig2.savefig(output_prefix + 'metabolite_validation.svg')
    
        fig3 = plt.figure(figsize=(6, 6))
        fig3.suptitle('Enzyme Concentrations')
        ax = fig3.add_subplot(1, 1, 1, xscale='log', yscale='log')
        model.ValidateEnzymeConcentrations(lnC_ECM, ax)
        fig3.savefig(output_prefix + 'enzyme_validation.svg')
    
    fig4 = plt.figure(figsize=(5, 5))
    ax = fig4.add_subplot(1, 1, 1)
    #model.PlotVolumesPie(lnC_ECM, ax)
    vols, labels, colors = model._GetVolumeDataForPlotting(lnC_ECM)
    ax.pie(vols, labels=labels, colors=colors)
    ax.set_title('total weight = %.2g [g/L]' % sum(vols))
    fig4.savefig(output_prefix + 'pie.svg')
