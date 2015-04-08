# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:38:13 2015

@author: noore
"""

from SBtabTools import oneOrMany
from SBtab import SBtabTable
from tablibIO import loadTSV
from python.component_contribution import ComponentContribution
from python.kegg_reaction import KeggReaction

class SBtabDict(dict):
    
    def __init__(self, sbtab_list):
        """
            Arguments:
                sbtab_list - a list of SBtabTable objects
        """
        for m in sbtab_list:
            self[m.table_name] = m

        self.reactions = self.GetColumnFromTable('Reaction', 'SumFormula')
        self.met2kegg = self.GetDictFromTable('Compound', 'Compound', 'Compound:Identifiers:kegg.compound')

    @staticmethod
    def FromSBtab(fpath):
        spreadsheet_file = loadTSV(fpath, False)
        m = oneOrMany(spreadsheet_file)
        sbtab_list = [SBtabTable(dset, fpath) for dset in m]
        return SBtabDict(sbtab_list)
        
    def GetColumnFromTable(self, table_name, column_name):
        column_index = self[table_name].columns_dict['!' + column_name]
        rows = self[table_name].getRows()
        return [r[column_index] for r in rows]
        
    def GetDictFromTable(self, table_name, key_column_name, value_column_name):
        key_index = self[table_name].columns_dict['!' + key_column_name]
        value_index = self[table_name].columns_dict['!' + value_column_name]
        rows = self[table_name].getRows()
        return {r[key_index] : r[value_index] for r in rows}