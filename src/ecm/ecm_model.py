# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:38:13 2015

@author: noore
"""

from SBtabTools import oneOrMany
from SBtab import SBtabTable
from tablibIO import loadTSV
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution import ComponentContribution
import re
import numpy as np
from component_contribution.thermodynamic_constants import default_RT as RT

class ECMmodel(object):
    
    def __init__(self, sbtab_fpath):
        self._sbtab_dict = SBtabDict.FromSBtab(sbtab_fpath)
        self.kegg_model = ECMmodel.GenerateKeggModel(self._sbtab_dict)
        
        self.kegg_model.check_S_balance()
        cc = ComponentContribution.init()
        self.kegg_model.add_thermo(cc)
        
        dG0_prime, sqrt_Sigma = self.kegg_model.get_transformed_dG0(pH=7.5, I=0.1, T=298.15)
        self.rid2dG0_cc = dict(zip(self.kegg_model.rids, dG0_prime.flat))
        
        rid2keq, rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2MM = \
            ECMmodel.ReadKineticParameters(self._sbtab_dict)
        
        self.rid2dG0_rate_constant_table = {rid: -RT*np.log(keq) for (rid, keq) in rid2keq.iteritems()}
        
        tmp = self._sbtab_dict.GetDictFromTable('GibbsEnergyOfReaction', 'Reaction', 'dG0')
        self.rid2dG0_gibbs_energy_table = {rid: float(dG0) for (rid, dG0) in tmp.iteritems()}

    @staticmethod
    def GenerateKeggModel(sbtab_dict):
        met2kegg = sbtab_dict.GetDictFromTable('Compound', 'Compound',
            'Compound:Identifiers:kegg.compound')
        
        reaction_names = sbtab_dict.GetColumnFromTable('Reaction', 'Reaction')
        reaction_formulas = sbtab_dict.GetColumnFromTable('Reaction', 'SumFormula')
        sparse_reactions = map(SBtabDict.ParseReaction, reaction_formulas)
        
        map_met2kegg = lambda spr : {met2kegg.get(k): v for (k,v) in spr.iteritems()}
        sparse_reactions_kegg = map(map_met2kegg, sparse_reactions)
        kegg_reactions = [KeggReaction(s, rid=rid) for (s, rid) 
                          in zip(sparse_reactions_kegg, reaction_names)]

        model = KeggModel.from_kegg_reactions(kegg_reactions, has_reaction_ids=True)
        return model

    @staticmethod
    def ReadKineticParameters(sbtab_dict, table_name='RateConstant'):
        cols = ['QuantityType',
                'Value',
                'Compound:Identifiers:kegg.compound',
                'Reaction',
                'Unit']

        rid2keq = {}       # equilibrium constants
        rid2crc_gmean = {} # catalytic rate constant geomertic mean
        rid2crc_fwd = {}   # catalytic rate constant forward
        rid2crc_rev = {}   # catalytic rate constant reverse
        crctype2dict = {'catalytic rate constant geometric mean': rid2crc_gmean,
                        'substrate catalytic rate constant': rid2crc_fwd,
                        'product catalytic rate constant': rid2crc_rev}
        
        rid_cid2MM = {}   # Michaelis-Menten constants
        
        for i, row in enumerate(sbtab_dict.GetColumnsFromTable(table_name, cols)):
            try:
                typ, val, cid, rid, unit = row
                val = float(val)
                
                if typ in crctype2dict:
                    if unit != '1/s':
                        raise AssertionError('Catalytic rate constants must be '
                                             'in units of 1/s, not %s' % unit)
                    crctype2dict[typ][rid] = val
                elif typ == 'equilibrium constant':
                    rid2keq[rid] = val
                elif typ == 'Michaelis constant':
                    if unit == 'mM':
                        rid_cid2MM[rid, cid] = val * 1e-3
                    elif unit == 'M':
                        rid_cid2MM[rid, cid] = val
                    else:
                        raise AssertionError('Michaelis constants must be in M or mM')
                else:
                    raise AssertionError('unrecognized Rate Constant Type: ' + typ)
            except AssertionError as e:
                raise ValueError('Syntax error in SBtab table %s, row %d - %s' %
                                 (table_name, i, str(e)))
                
        return rid2keq, rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2MM
        
class SBtabDict(dict):
    
    def __init__(self, sbtab_list):
        """
            Arguments:
                sbtab_list - a list of SBtabTable objects
        """
        for m in sbtab_list:
            self[m.table_name] = m

    @staticmethod
    def FromSBtab(fpath):
        spreadsheet_file = loadTSV(fpath, False)
        m = oneOrMany(spreadsheet_file)
        sbtab_list = [SBtabTable(dset, fpath) for dset in m]
        return SBtabDict(sbtab_list)

    def GetColumnFromTable(self, table_name, column_name):
        """
            Returns:
                a list of the values in the column called 'column_name'
                in the table 'table_name'
        """
        column_index = self[table_name].columns_dict['!' + column_name]
        rows = self[table_name].getRows()
        return [r[column_index] for r in rows]

    def GetColumnsFromTable(self, table_name, column_names):
        """
            Arguments:
                table_name   - the name of the table in the SBtab file (without '!!')
                column_names - a list of column names from which to get the data (without '!')
                
            Returns:
                a list of lists containing the values corresponding to the
                columns in 'column_names' in the table 'table_name'
        """
        idxs = [self[table_name].columns_dict['!' + c] for c in column_names]
        return [map(r.__getitem__, idxs) for r in self[table_name].getRows()]
        
    def GetDictFromTable(self, table_name, key_column_name, value_column_name):
        column_names = [key_column_name, value_column_name]        
        return dict(self.GetColumnsFromTable(table_name, column_names))

    @staticmethod
    def ParseReactionFormulaSide(s):
        """ 
            Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
            Ignores stoichiometry.
            
            Returns:
                The set of CIDs.
        """
        if s.strip() == "null":
            return {}
        
        compound_bag = {}
        for member in re.split('\s+\+\s+', s):
            tokens = member.split(None, 1)
            if len(tokens) == 0:
                continue
            if len(tokens) == 1:
                amount = 1
                key = member
            else:
                try:
                    amount = float(tokens[0])
                except ValueError:
                    raise Exception(
                        "Non-specific reaction: %s" % s)
                key = tokens[1]
                
            try:
                compound_bag[key] = compound_bag.get(key, 0) + amount
            except ValueError:
                raise Exception(
                    "Non-specific reaction: %s" % s)
        
        return compound_bag

    @staticmethod
    def ParseReaction(formula, arrow='<=>'):
        """ 
            Parses a two-sided formula such as: 2 FBP => DHAP + GAP
            
            Return:
                The set of substrates, products and the direction of the reaction
        """
        tokens = formula.split(arrow)
        if len(tokens) < 2:
            raise Exception('Reaction does not contain the arrow sign (%s): %s'
                            % (arrow, formula))
        if len(tokens) > 2:
            raise Exception('Reaction contains more than one arrow sign (%s): %s'
                            % (arrow, formula))
        
        left = tokens[0].strip()
        right = tokens[1].strip()
        
        sparse_reaction = {}
        for cid, count in SBtabDict.ParseReactionFormulaSide(left).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 

        for cid, count in SBtabDict.ParseReactionFormulaSide(right).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 

        return sparse_reaction
    
