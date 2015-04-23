# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:38:13 2015

@author: noore
"""

from sbtab_dict import SBtabDict
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution import ComponentContribution
import numpy as np
from cost_function import EnzymeCostFunction
from scipy.io import savemat
import colors
from util import ParseReaction, PlotCorrelation
import logging
from errors import ThermodynamicallyInfeasibleError

CELL_VOL_PER_DW = 2.7e-3 # L/gCDW [Winkler and Wilson, 1966, http://www.jbc.org/content/241/10/2200.full.pdf+html]

class ECMmodel(object):
    
    def __init__(self, sbtab_fpath):
        self._sbtab_dict = SBtabDict.FromSBtab(sbtab_fpath)
        self.kegg2met = self._sbtab_dict.GetDictFromTable('Compound', 
            'Compound:Identifiers:kegg.compound', 'NameForPlots')
        self.kegg2rxn = self._sbtab_dict.GetDictFromTable('Reaction', 
            'Reaction', 'NameForPlots')
        self.kegg_model = ECMmodel.GenerateKeggModel(self._sbtab_dict)

        # a dictionary indicating which compound is external or not
        self.cid2external = self._sbtab_dict.GetDictFromTable(
            'Compound', 'Compound:Identifiers:kegg.compound', 'External',
            value_mapping=bool)

        # read the lower and upper bounds and convert them to M
        bound_units = self._sbtab_dict.GetTableAttribute('ConcentrationConstraint', 'Unit')
        bound_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(bound_units)
        self.cid2min_bound = self._sbtab_dict.GetDictFromTable(
            'ConcentrationConstraint', 'Compound:Identifiers:kegg.compound', 'Concentration:Min',
            value_mapping=bound_mapping)
        self.cid2max_bound = self._sbtab_dict.GetDictFromTable(
            'ConcentrationConstraint', 'Compound:Identifiers:kegg.compound', 'Concentration:Max',
            value_mapping=bound_mapping)

        self.kegg_model.check_S_balance()
        cc = ComponentContribution.init()
        self.kegg_model.add_thermo(cc)
        self._CalcGibbsEneriges()
        
        rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2KMM = \
            ECMmodel._ReadKineticParameters(self._sbtab_dict)
        
        # read flux values and convert them to M/s
        flux_units = self._sbtab_dict.GetTableAttribute('Flux', 'Unit')
        flux_mapping = ECMmodel._MappingToCanonicalFluxUnits(flux_units)
        self.rid2flux = self._sbtab_dict.GetDictFromTable(
            'Flux', 'Reaction', 'Flux', value_mapping=flux_mapping)
        
        S = self.kegg_model.S
        flux = np.matrix(map(self.rid2flux.get, self.kegg_model.rids)).T
        kcat = np.matrix(map(rid2crc_gmean.get, self.kegg_model.rids)).T
        dG0 = np.matrix(map(self.rid2dG0.get, self.kegg_model.rids)).T
        KMM = ECMmodel._GenerateKMM(self.kegg_model.cids,
                                    self.kegg_model.rids, rid_cid2KMM)
        c_bounds = np.array(zip(map(self.cid2min_bound.get, self.kegg_model.cids),
                                map(self.cid2max_bound.get, self.kegg_model.cids)))
        lnC_bounds = np.log(c_bounds) # assume bounds are in M
        self.ecf = EnzymeCostFunction(S, flux=flux, kcat=kcat, dG0=dG0, KMM=KMM,
                                      lnC_bounds=lnC_bounds, ecf_version='ECF3')
    
    def WriteMatFile(self, file_name):
        mdict = self.ecf.Serialize()
        mdict['cids'] = self.kegg_model.cids
        mdict['rids'] = self.kegg_model.rids
        savemat(file_name, mdict, format='5')
    
    @staticmethod
    def GenerateKeggModel(sbtab_dict):
        met2kegg = sbtab_dict.GetDictFromTable('Compound', 'Compound',
            'Compound:Identifiers:kegg.compound')
        
        reaction_names = sbtab_dict.GetColumnFromTable('Reaction', 'Reaction')
        reaction_formulas = sbtab_dict.GetColumnFromTable('Reaction', 'SumFormula')
        sparse_reactions = map(ParseReaction, reaction_formulas)
        
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

        rid2crc_gmean = {} # catalytic rate constant geomertic mean [1/s]
        rid2crc_fwd = {}   # catalytic rate constant forward [1/s]
        rid2crc_rev = {}   # catalytic rate constant reverse [1/s]
        crctype2dict = {'catalytic rate constant geometric mean': rid2crc_gmean,
                        'substrate catalytic rate constant': rid2crc_fwd,
                        'product catalytic rate constant': rid2crc_rev}
        
        rid_cid2KMM = {}   # Michaelis-Menten constants [M]
        
        for i, row in enumerate(sbtab_dict.GetColumnsFromTable(table_name, cols)):
            try:
                typ, val, cid, rid, unit = row
                val = float(val)
                
                if typ in crctype2dict:
                    if unit != '1/s':
                        raise AssertionError('Catalytic rate constants must be '
                                             'in units of 1/s, not %s' % unit)
                    crctype2dict[typ][rid] = val
                elif typ == 'Michaelis constant':
                    value_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(unit)
                    rid_cid2KMM[rid, cid] = value_mapping(val)
                else:
                    raise AssertionError('unrecognized Rate Constant Type: ' + typ)
            except AssertionError as e:
                raise ValueError('Syntax error in SBtab table %s, row %d - %s' %
                                 (table_name, i, str(e)))
                
        return rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2KMM

    def _CalcGibbsEneriges(self):
        dG0_prime, sqrt_Sigma = self.kegg_model.get_transformed_dG0(pH=7.5, I=0.1, T=298.15)
        self.rid2dG0 = dict(zip(self.kegg_model.rids, dG0_prime.flat))

        # legacy code - read the dG0 from the SBtab itself rather than calculating
        # it using CC
        #
        #self.rid2dG0 = self._sbtab_dict.GetDictFromTable(
        #    'GibbsEnergyOfReaction', 'Reaction', 'dG0', value_mapping=float)
    
    @staticmethod
    def _GenerateKMM(cids, rids, rid_cid2KMM):
        KMM = np.ones((len(cids), len(rids)))
        for i, cid in enumerate(cids):
            for j, rid in enumerate(rids):
                KMM[i, j] = rid_cid2KMM.get((rid,cid), 1)
        return KMM
    
    def MDF(self):
        mdf, params = self.ecf.MDF()
        if np.isnan(mdf) or mdf < 0.0:
            logging.error('Negative MDF value: %.1f' % mdf)
            logging.error('The reactions with shadow prices are:')
            shadow_prices = params['reaction prices']
            for rid, sp in zip(self.kegg_model.rids, shadow_prices.flat):
                if sp:
                    logging.error('\t%s : %g' % (rid, sp))
            raise ThermodynamicallyInfeasibleError()
        return params['ln concentrations']
        
    def ECM(self, lnC0=None):
        return self.ecf.ECM(lnC0 or self.MDF())
        
    def ECF(self, lnC):
        return self.ecf.ECF(lnC)
        
    @staticmethod
    def _MappingToCanonicalConcentrationUnits(unit):
        """
            Assuming the canonical units for concentration are Molar
            
            Returns:
                A function that converts a single number or string to the 
                canonical units
        """
        if unit == 'M':
            return lambda x: float(x)
        if unit == 'mM':
            return lambda x: float(x)*1e-3
        if unit == 'uM':
            return lambda x: float(x)*1e-6
        if unit == 'nM':
            return lambda x: float(x)*1e-9
        
        raise ValueError('Cannot convert these units to M: ' + unit)
    
    @staticmethod
    def _MappingToCanonicalFluxUnits(unit):
        """
            Assuming the canonical units for concentration are [M/s]
            
            Returns:
                A function that converts a single number or string to the 
                canonical units
                
            Note: CELL_VOL_PER_DW is given in [L/gCDW]
        """
        if unit == 'M/s':
            return lambda x: float(x)
        if unit == 'mM/s':
            return lambda x: float(x)*1e-3
        if unit == 'mmol/gCDW/h':
            return lambda x: float(x) / (CELL_VOL_PER_DW * 3600 * 1e3)
        if unit == 'mol/gCDW/h':
            return lambda x: float(x) / (CELL_VOL_PER_DW * 3600)
        if unit == 'umol/gCDW/min':
            return lambda x: float(x) / (CELL_VOL_PER_DW * 60 * 1e6)
        if unit == 'mmol/gCDW/min':
            return lambda x: float(x) / (CELL_VOL_PER_DW * 60 * 1e3)
        
        raise ValueError('Cannot convert these units to M/s: ' + unit)

    def _GetMeasuredMetaboliteConcentrations(self):
        unit = self._sbtab_dict.GetTableAttribute('Concentration', 'Unit')
        value_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(unit)
        
        # assume concentrations are in mM
        return self._sbtab_dict.GetDictFromTable(
            'Concentration', 'Compound:Identifiers:kegg.compound',
            'Concentration', value_mapping=value_mapping)
        
    def _GetMeasuredEnzymeConcentrations(self):
        unit = self._sbtab_dict.GetTableAttribute('EnzymeConcentration', 'Unit')
        value_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(unit)

        return self._sbtab_dict.GetDictFromTable(
            'EnzymeConcentration', 'Reaction',
            'EnzymeConcentration', value_mapping=value_mapping)

    def PlotEnzymeCosts(self, lnC, ax, top_level=3):
        """
            A bar plot in log-scale showing the partitioning of cost between
            the levels of kinetic costs:
            1 - capacity
            2 - thermodynamics
            3 - saturation
            4 - allosteric
        """
        assert top_level in range(1, 5)
        
        ecf_mat = self.ecf.GetEnzymeCostPartitions(lnC)
        datamat = np.log(ecf_mat)
        base = min(datamat[np.isfinite(datamat)].flat) - 1
        bottoms = np.hstack([np.ones((datamat.shape[0], 1)) * base,
                             np.cumsum(datamat, 1)])
        bottoms = np.exp(bottoms)
        steps = np.diff(bottoms)
        
        labels = ['capacity', 'thermodynamic', 'saturation', 'allosteric']
        labels = labels[0:top_level]

        ind = np.arange(ecf_mat.shape[0])    # the x locations for the groups
        width = 0.7
        cmap = colors.ColorMap(labels, saturation=0.5, value=0.8)
        ax.set_yscale('log')
        for i, label in enumerate(labels):
            ax.bar(ind, steps[:, i], width,
                   bottom=bottoms[:, i], color=cmap[label],
                   alpha=1.0)
        ax.set_xticks(ind + width/2)
        xticks = map(self.kegg2rxn.get, self.kegg_model.rids)
        ax.set_xticklabels(xticks, size='medium', rotation=45)
        ax.legend(labels, loc='best', framealpha=0.2)
        #ax.set_xlabel('reaction')
        ax.set_ylabel('enzyme cost [M]')
        ax.set_ylim(ymin=base)
        total = np.prod(ecf_mat, 1).sum()
        ax.set_title(r'Total enzyme cost = %.2f $\times$ $10^{-3}$ [M]' % (total*1e3))
    
    def ValidateMetaboliteConcentrations(self, lnC, ax):
        pred_conc = np.exp(lnC)

        meas_met2conc = self._GetMeasuredMetaboliteConcentrations()
        meas_conc = np.matrix(map(meas_met2conc.get, self.kegg_model.cids)).T
        
        mask =  (meas_conc > 0) & (pred_conc > 0)   # remove NaNs and zeros
        mask &= np.diff(self.ecf.lnC_bounds) > 1e-9 # remove compounds with fixed concentrations

        labels = map(self.kegg2met.get, self.kegg_model.cids)
        PlotCorrelation(ax, meas_conc, pred_conc, labels, mask)
        ax.set_xlabel('measured [M]')
        ax.set_ylabel('predicted [M]')
        
    def ValidateEnzymeConcentrations(self, lnC, ax):
        pred_conc = self.ecf.ECF(lnC)

        meas_enz2conc = self._GetMeasuredEnzymeConcentrations()
        meas_conc = np.matrix(map(meas_enz2conc.get, self.kegg_model.rids)).T
        
        mask = (pred_conc > 0) & (meas_conc > 0)

        labels = map(self.kegg2rxn.get, self.kegg_model.rids)
        PlotCorrelation(ax, meas_conc, pred_conc, labels, mask)

        ax.set_xlabel('measured [M]')
        ax.set_ylabel('predicted [M]')

    def WriteHtmlTables(self, lnC, html):
        met_conc = dict(zip(self.kegg_model.cids, np.exp(lnC).flat))
        cid2lower_bound = dict(zip(self.kegg_model.cids, np.exp(self.ecf.lnC_bounds[:, 0].flat)))
        cid2upper_bound = dict(zip(self.kegg_model.cids, np.exp(self.ecf.lnC_bounds[:, 1].flat)))
        enz_conc = dict(zip(self.kegg_model.rids, self.ecf.ECF(lnC).flat))
        driving_forces = dict(zip(self.kegg_model.rids, self.ecf._DrivingForce(lnC).flat))
        
        headers = ['Reaction', 'KEGG ID', 'flux [mM/s]', 'enzyme conc. [uM]',
                   'Driving force [kJ/mol]']
        values = [(self.kegg2rxn[r], r, self.rid2flux[r]*1e3, enz_conc[r]*1e6,
                   driving_forces[r])
                  for r in self.kegg_model.rids]
        rowdicst = [dict(zip(headers, v)) for v in values]
        html.write_table(rowdicst, headers=headers, decimal=3)
        
        headers = ['Compound', 'KEGG ID', 'Concentration [M]',
                   'Lower bound [M]', 'Upper bound [M]']
        values = [(self.kegg2met[cid], cid, '%.2e' % met_conc[cid],
                   '%.2e' % cid2lower_bound[cid], '%.2e' % cid2upper_bound[cid])
                  for cid in self.kegg_model.cids]
        rowdicst = [dict(zip(headers, v)) for v in values]
        html.write_table(rowdicst, headers=headers)

