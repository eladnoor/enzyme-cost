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

def MakeParser():
    parser = argparse.ArgumentParser(
        description=('Run Enzyme Cost Minimization (ECM) on a model '
                     'given as a single SBtab file'))
    parser.add_argument('sbtab', type=str, help='SBtab input filename')
    parser.add_argument('respath', type=str, help='path for writing result files')
    
    parser.add_argument('--level', type=int, help='Enzyme Cost Function level: 1, 2, [3], or 4',
                        default=3)
    parser.add_argument('--dgsource', type=str,
                        help="Source for dG0s: [keq_table], dG0r_table, or component_contribution",
                        default='keq_table')
    parser.add_argument('--kcatsource', type=str,
                        help="Source for kcats: fwd or [gmean]",
                        default='gmean')
    parser.add_argument('--denominator', type=str,
                        help="Rate law denominator: S, SP, 1S, 1SP, or [CM]",
                        default='CM')
    parser.add_argument('--regularization', type=str,
                        help="Regularization method: none, [volume], or quadratic",
                        default='volume')
    return parser

###############################################################################
parser = MakeParser()
args = parser.parse_args()
logging.getLogger().setLevel(logging.WARNING)

SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
DATA_DIR = os.path.join(os.path.split(SCRIPT_DIR)[0], 'data')
RESULT_DIR = os.path.join(os.path.split(SCRIPT_DIR)[0], 'res')

logging.info('Reading SBtab files')
modeldata_sbtabs = SBtabDict.FromSBtab(args.sbtab)

logging.info('Creating an ECM model using the data')

ecf_params = {
    'version'       : args.level,
    'dG0_source'    : args.dgsource,
    'kcat_source'   : args.kcatsource,
    'denominator'   : args.denominator,
    'regularization': args.regularization
    }

model = ECMmodel(modeldata_sbtabs, ecf_params=ecf_params)

logging.info('Solving MDF problem')
lnC_MDF = model.MDF()
logging.info('Solving ECM problem')
lnC_ECM = model.ECM(n_iter=5)

res_sbtab = model.ToSBtab(lnC_ECM,
                          args.respath,
                          document_name='E. coli central carbon metabolism - ECM result')
