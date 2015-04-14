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
from cost_function import EnzymeCostFunction

class ECMmodel(object):
    
    def __init__(self, sbtab_fpath):
        self._sbtab_dict = SBtabDict.FromSBtab(sbtab_fpath)
        self.kegg_model = ECMmodel.GenerateKeggModel(self._sbtab_dict)

        # a dictionary indicating which compound is external or not
        self.cid2external = self._sbtab_dict.GetDictFromTable(
            'Compound', 'Compound:Identifiers:kegg.compound', 'External',
            value_mapping=bool)

        # read the lower and upper bounds. assume they are given in mM
        self.cid2min_bound = self._sbtab_dict.GetDictFromTable(
            'ConcentrationConstraint', 'Compound:Identifiers:kegg.compound', 'Concentration:Min',
            value_mapping=float)
        self.cid2max_bound = self._sbtab_dict.GetDictFromTable(
            'ConcentrationConstraint', 'Compound:Identifiers:kegg.compound', 'Concentration:Max',
            value_mapping=float)

        self.kegg_model.check_S_balance()
        cc = ComponentContribution.init()
        self.kegg_model.add_thermo(cc)
        self._CalcGibbsEneriges()
        
        rid2keq, rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2KMM = \
            ECMmodel._ReadKineticParameters(self._sbtab_dict)
        
        # load the standard Gibbs energies from the SBtab file (for legacy reasons)
        # basically, we don't need them because we can get the same data
        # from component-contribution directly
        
        #rid2dG0_rate_constant_table = {rid: -RT*np.log(keq) for (rid, keq) in rid2keq.iteritems()}
        #
        #rid2dG0_gibbs_energy_table = self._sbtab_dict.GetDictFromTable(
        #    'GibbsEnergyOfReaction', 'Reaction', 'dG0', value_mapping=float)
        
        rid2flux = self._sbtab_dict.GetDictFromTable(
            'Flux', 'Reaction', 'Flux', value_mapping=float)
        
        S = self.kegg_model.S
        flux = np.matrix(map(rid2flux.get, self.kegg_model.rids)).T
        kcat = np.matrix(map(rid2crc_gmean.get, self.kegg_model.rids)).T
        dG0 = np.matrix(map(self.rid2dG0.get, self.kegg_model.rids)).T
        KMM = ECMmodel._GenerateKMM(self.kegg_model.cids,
                                    self.kegg_model.rids, rid_cid2KMM)
        c_bounds = np.array(zip(map(self.cid2min_bound.get, self.kegg_model.cids),
                                map(self.cid2max_bound.get, self.kegg_model.cids)))
        lnC_bounds = np.log(c_bounds * 1e-3) # convert from mM to M
        self.ecf = EnzymeCostFunction(S, flux=flux, kcat=kcat, dG0=dG0, KMM=KMM,
                                      lnC_bounds=lnC_bounds, ecf_version='ECF2')
    
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
    def _ReadKineticParameters(sbtab_dict, table_name='RateConstant'):
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
        
        rid_cid2KMM = {}   # Michaelis-Menten constants
        
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
                        rid_cid2KMM[rid, cid] = val * 1e-3
                    elif unit == 'M':
                        rid_cid2KMM[rid, cid] = val
                    else:
                        raise AssertionError('Michaelis constants must be in M or mM')
                else:
                    raise AssertionError('unrecognized Rate Constant Type: ' + typ)
            except AssertionError as e:
                raise ValueError('Syntax error in SBtab table %s, row %d - %s' %
                                 (table_name, i, str(e)))
                
        return rid2keq, rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2KMM

    def _CalcGibbsEneriges(self, mode='CC'):
        if mode == 'CC':
            dG0_prime, sqrt_Sigma = self.kegg_model.get_transformed_dG0(pH=7.5, I=0.1, T=298.15)
            self.rid2dG0 = dict(zip(self.kegg_model.rids, dG0_prime.flat))
        elif mode == 'KEQ':
            # loading the standard Gibbs energies from the SBtab file (for legacy reasons)
            # basically, we don't need them because we can get the same data
            # from component-contribution directly

            rid2keq, _, _, _, _ = ECMmodel.ReadKineticParameters(self._sbtab_dict)
            self.rid2dG0 = {rid: -RT*np.log(keq) for (rid, keq) in rid2keq.iteritems()}
        elif mode == 'DG0':
            self.rid2dG0 = self._sbtab_dict.GetDictFromTable(
                'GibbsEnergyOfReaction', 'Reaction', 'dG0', value_mapping=float)
    
    @staticmethod
    def _GenerateKMM(cids, rids, rid_cid2KMM):
        KMM = np.ones((len(cids), len(rids)))
        for i, cid in enumerate(cids):
            for j, rid in enumerate(rids):
                KMM[i, j] = rid_cid2KMM.get((rid,cid), 1)
        return KMM
    
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
        
    def GetDictFromTable(self, table_name, key_column_name, value_column_name,
                         value_mapping=None):
        column_names = [key_column_name, value_column_name]   
        keys, vals = zip(*self.GetColumnsFromTable(table_name, column_names))
        return dict(zip(keys, map(value_mapping, vals)))

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
    
