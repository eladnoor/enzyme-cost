# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:38:13 2015

@author: noore
"""

import logging
import numpy as np
from collections import defaultdict

from ecm.kegg_model import KeggModel, KeggReaction, ParseReaction
from ecm.cost_function import EnzymeCostFunction
from scipy.io import savemat
from ecm.colors import ColorMap
from ecm.util import RT, CELL_VOL_PER_DW, ECF_DEFAULTS, str2bool, PlotCorrelation
from ecm.errors import ThermodynamicallyInfeasibleError

import sys, os
import pandas as pd

class ECMmodel(object):

    DATAFRAME_NAMES = ['Compound', 'Reaction', 'ConcentrationConstraint',
                       'Parameter', 'Flux']
    
    def __init__(self, df_dict, bound_unit='mM',
                 flux_unit='mM/s', validate_sbtab=None,
                 ecf_params=None):
        
        for n in ECMmodel.DATAFRAME_NAMES:
            if n not in df_dict.keys():
                raise KeyError('A DataFrame named %s is missing from the model' % n)
        
        self.ecf_params = dict(ECF_DEFAULTS)
        if ecf_params is not None:
            self.ecf_params.update(ecf_params)
        
        self._validate_sbtab = validate_sbtab
        
        if 'NameForPlots' in df_dict['Compound']:
            self.kegg2met = df_dict['Compound'].set_index('Identifiers:kegg.compound')['NameForPlots']
        elif 'Name' in df_dict['Compound']:
            self.kegg2met = df_dict['Compound'].set_index('Identifiers:kegg.compound')['Name']
        else:
            raise KeyError('Column "Name" or "NameForPlots" must be given for the Compounds table')
            
        reaction_df = df_dict['Reaction']
        if 'NameForPlots' in reaction_df:
            self.kegg2rxn = reaction_df.set_index('ID')['NameForPlots']
        elif 'Name' in reaction_df:
            self.kegg2rxn = reaction_df.set_index('ID')['Name']
        else:
            raise KeyError('Column "Name" or "NameForPlots" must be given for the Reactions table')

        self.kegg_model = ECMmodel.GenerateKeggModel(df_dict['Compound'],
                                                     df_dict['Reaction'])

        # a dictionary indicating which compound is external or not
        if 'External' in df_dict['Compound']:
            self.cid2external = df_dict['Compound']['External'].apply(str2bool)
        elif 'IsConstant' in df_dict['Compound']:
            self.cid2external = df_dict['Compound']['IsConstant'].apply(str2bool)
        else:
            raise KeyError('Column "External" or "IsConstant" must be given for the Compounds table')

        self.cid2min_bound, self.cid2max_bound = \
            ECMmodel.ReadConcentrationBounds(bound_unit, df_dict['ConcentrationConstraint'])

        try:
            from component_contribution.kegg_model import KeggModel as CCKeggModel
            self.cc_model = CCKeggModel(self.kegg_model.S.copy(),
                                        list(self.kegg_model.cids),
                                        list(self.kegg_model.rids))
            self.cc_model.check_S_balance(fix_water=True)
        except ImportError:
            self.cc_model = None

        rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2KMM, rid2keq, rid2mw, cid2mw = \
            ECMmodel.ReadParameters(df_dict['Parameter'])

        if self.ecf_params['dG0_source'] == 'keq_table':
            self._CalcGibbsEnergiesFromKeq(rid2keq)
        elif self.ecf_params['dG0_source'] == 'component_contribution':
            self._CalcGibbsEnerigesFromComponentContribution()
        elif self.ecf_params['dG0_source'] == 'dG0r_table':
            self._CalcGibbsEnerigesFromdG0ReactionTable()
        else:
            raise ValueError('unrecognized dG0 source: ' +
                             self.ecf_params['dG0_source'])

        # read flux values and convert them to M/s
        flux_mapping = ECMmodel._MappingToCanonicalFluxUnits(flux_unit)
        self.rid2flux = df_dict['Flux'].set_index('Reaction')
        
        flux = df_dict['Flux'].set_index('Reaction').loc[self.kegg_model.rids, 'Value']
        flux = np.matrix(flux.apply(flux_mapping).values).T
        
        # we only need to define get kcat in the direction of the flux
        # if we use the 'gmean' option, that means we assume we only know
        # the geometric mean of the kcat, and we distribute it between
        # kcat_fwd and kcat_bwd according to the Haldane relationship
        # if we use the 'fwd' option, we just take the kcat in the
        # direction of flux (as is) and that would mean that our
        # thermodynamic rate law would be equivalent to calculating the
        # reverse kcat using the Haldane relationship
        if self.ecf_params['kcat_source'] == 'gmean':
            kcat = np.matrix(list(map(rid2crc_gmean.get, self.kegg_model.rids))).T
        elif self.ecf_params['kcat_source'] == 'fwd':
            # get the relevant kcat (fwd/rev) depending on the direction of flux
            kcat = []
            for r, rid in enumerate(self.kegg_model.rids):
                if flux[r] >= 0:
                    kcat.append(rid2crc_fwd[rid])
                else:
                    kcat.append(rid2crc_rev[rid])
            kcat = np.matrix(kcat).T
        else:
            raise ValueError('unrecognized kcat source: ' +
                             self.ecf_params['kcat_source'])

        dG0 = np.matrix(list(map(self.rid2dG0.get, self.kegg_model.rids))).T
        KMM = ECMmodel._GenerateKMM(self.kegg_model.cids,
                                    self.kegg_model.rids, rid_cid2KMM)
        c_bounds = np.vstack([self.cid2min_bound[self.kegg_model.cids].values,
                              self.cid2max_bound[self.kegg_model.cids].values]).T
        lnC_bounds = np.log(c_bounds) # assume bounds are in M

        # we need all fluxes to be positive, so for every negative flux,
        # we multiply it and the corresponding column in S by (-1)
        dir_mat = np.matrix(np.diag(np.sign(flux + 1e-10).flat))
        flux = dir_mat * flux
        S = self.kegg_model.S * dir_mat
        dG0 = dir_mat * dG0

        mw_enz = np.matrix(list(map(rid2mw.get, self.kegg_model.rids))).T
        mw_met = np.matrix(list(map(cid2mw.get, self.kegg_model.cids))).T

        self.ecf = EnzymeCostFunction(S, flux=flux, kcat=kcat,
                                      dG0=dG0, KMM=KMM,
                                      mw_enz=mw_enz, mw_met=mw_met,
                                      lnC_bounds=lnC_bounds,
                                      params=self.ecf_params)
    
    def from_sbtabs(model_sbtab, 
                    validate_sbtab=None, ecf_params=None):
        
        df_dict = {k: v.toDataFrame() for k, v in model_sbtab.items()}
        
        bound_unit = model_sbtab.GetTableAttribute('ConcentrationConstraint', 'Unit')
        flux_unit = model_sbtab.GetTableAttribute('Flux', 'Unit')
        return ECMmodel(df_dict,
                        bound_unit=bound_unit,
                        flux_unit=flux_unit,
                        validate_sbtab=validate_sbtab,
                        ecf_params=ecf_params)

    def WriteMatFile(self, file_name):
        mdict = self.ecf.Serialize()
        mdict['cids'] = self.kegg_model.cids
        mdict['rids'] = self.kegg_model.rids
        with open(file_name, 'wb') as fp:
            savemat(fp, mdict, format='5')

    @staticmethod
    def GenerateKeggModel(compound_df, reaction_df):
        # first, parse the reaction formulas and create a dictionary
        # for each one (i.e. "sparse reaction")
        sparse_reactions = list(map(ParseReaction,
                                    reaction_df['ReactionFormula']))
        
        # use the compound table to convert the IDs to KEGG IDs
        met2kegg = compound_df.set_index('ID')['Identifiers:kegg.compound']
        map_met2kegg = lambda spr : {met2kegg[k]: v for (k,v) in spr.items()}
        try:
            sparse_reactions_kegg = list(map(map_met2kegg, sparse_reactions))
        except KeyError as e:
            raise KeyError('One of the metabolites in the reaction formulas '
                            'is not defined in the Compound table: ' + str(e))

        kegg_reactions = [KeggReaction(s, rid=rid) for (s, rid)
                          in zip(sparse_reactions_kegg, reaction_df['ID'])]

        model = KeggModel.from_kegg_reactions(kegg_reactions,
                                              has_reaction_ids=True)
        return model

    @staticmethod
    def ReadParameters(parameter_df):
        cols = ['QuantityType',
                'Value',
                'Compound:Identifiers:kegg.compound',
                'Reaction',
                'Unit']

        rid2mw = defaultdict(float)
        cid2mw = defaultdict(float)
        rid2keq = {}
        rid2crc_gmean = {} # catalytic rate constant geomertic mean
        rid2crc_fwd = {}   # catalytic rate constant forward
        rid2crc_rev = {}   # catalytic rate constant reverse
        crctype2dict = {'catalytic rate constant geometric mean': rid2crc_gmean,
                        'substrate catalytic rate constant': rid2crc_fwd,
                        'product catalytic rate constant': rid2crc_rev}

        rid_cid2KMM = {}   # Michaelis-Menten constants

        for i, row in parameter_df.iterrows():
            try:
                typ, val, cid, rid, unit = [row[c] for c in cols]
                val = float(val)

                if typ in crctype2dict:
                    if unit != '1/s':
                        raise AssertionError('Catalytic rate constants must be '
                                             'in units of 1/s, not %s' % unit)
                    crctype2dict[typ][rid] = val
                elif typ == 'Michaelis constant':
                    value_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(unit)
                    rid_cid2KMM[rid, cid] = value_mapping(val)
                elif typ == 'equilibrium constant':
                    assert unit == 'dimensionless'
                    rid2keq[rid] = val
                elif typ == 'protein molecular mass':
                    value_mapping = ECMmodel._MappingToCanonicalMolecularWeightUnits(unit)
                    rid2mw[rid] = value_mapping(val)
                elif typ == 'molecular mass':
                    value_mapping = ECMmodel._MappingToCanonicalMolecularWeightUnits(unit)
                    cid2mw[cid] = value_mapping(val)
                else:
                    raise AssertionError('unrecognized Rate Constant Type: ' + typ)
            except AssertionError:
                raise ValueError('Syntax error in Parameter table, row %d - %s' %
                                 (i, row))
        # make sure not to count water as contributing to the volume or
        # cost of a reaction
        return rid2crc_gmean, rid2crc_fwd, rid2crc_rev, rid_cid2KMM, rid2keq, rid2mw, cid2mw

    def _GibbsEnergyFromMilliMolarToMolar(self, dGm):
        """
            In the current SBtab file format, the 'standard Gibbs energies'
            are actually dGm (since mM is used as the reference concentration).
            We need to convert them back to M and that requires using the
            stoichiometric matrix (self.kegg_model.S).
        """

        # Assume that all concentrations are 1 mM
        mM_conc_v = 1e-3 * np.matrix(np.ones((self.kegg_model.S.shape[0], 1)))

        # subtract the effect of the concentrations to return back to dG0
        dG0 = dGm - self.kegg_model.S.T * RT * np.log(mM_conc_v)

        return dG0

    def _CalcGibbsEnergiesFromKeq(self, rid2keq):
        keq = np.matrix(list(map(rid2keq.get, self.kegg_model.rids))).T
        dGm_prime = - RT * np.log(keq)
        dG0_prime = self._GibbsEnergyFromMilliMolarToMolar(dGm_prime)
        self.rid2dG0 = dict(zip(self.kegg_model.rids, dG0_prime.flat))

    def _CalcGibbsEnerigesFromComponentContribution(self):
        try:
            from component_contribution.component_contribution_trainer import ComponentContribution
        except ImportError:
            raise Exception('You must install component_contribution in order '
                            'to calculate the reaction Gibbs energies: '
                            'https://github.com/eladnoor/component-contribution')
        # convert the kegg_model to a "real" KeggModel, i.e. one that
        # supports component-contributions
        cc = ComponentContribution.init()
        self.cc_model.add_thermo(cc)
        dG0_prime, dG0_std, sqrt_Sigma = self.cc_model.get_transformed_dG0(
            pH=7.5, I=0.1, T=298.15)
        self.rid2dG0 = dict(zip(self.kegg_model.rids, dG0_prime.flat))

    def _CalcGibbsEnerigesFromdG0ReactionTable(self):
        """
            it's better to take the values from the SBtab since they are processed
            by the parameter balancing funcion
        """
        dG0_units = self._model_sbtab.GetTableAttribute('GibbsEnergyOfReaction', 'Unit')
        value_mapping = ECMmodel._MappingToCanonicalEnergyUnits(dG0_units)
        rid2dGm = self._model_sbtab.GetDictFromTable(
            'GibbsEnergyOfReaction', 'Reaction', 'Value', value_mapping=value_mapping)
        dGm_prime = np.matrix(list(map(rid2dGm.get, self.kegg_model.rids), dtype=float)).T
        dG0_prime = self._GibbsEnergyFromMilliMolarToMolar(dGm_prime)
        self.rid2dG0 = dict(zip(self.kegg_model.rids, dG0_prime.flat))

    @staticmethod
    def ReadConcentrationBounds(bound_unit, bound_df):
        """
            read the lower and upper bounds and convert them to M.
            verify that the values make sense.
        """
        bound_df = bound_df.set_index('Compound:Identifiers:kegg.compound')
        bound_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(bound_unit)
        cid2min_bound = bound_df['Concentration:Min'].apply(bound_mapping)
        cid2max_bound = bound_df['Concentration:Max'].apply(bound_mapping)

        for cid in cid2min_bound.keys():
            assert cid2min_bound[cid] <= cid2max_bound[cid]
        return cid2min_bound, cid2max_bound

    @staticmethod
    def _GenerateKMM(cids, rids, rid_cid2KMM):
        KMM = np.ones((len(cids), len(rids)))
        for i, cid in enumerate(cids):
            for j, rid in enumerate(rids):
                KMM[i, j] = rid_cid2KMM.get((rid,cid), 1)
        return KMM

    def MDF(self):
        mdf, params = self.ecf.MDF() # note that MDF is given in units of RT

        if np.isnan(mdf) or mdf < 0.0:
            logging.error('Negative MDF value: %.1f RT (= %.1f kJ/mol)' %
                          (mdf, mdf * RT))
            logging.error('The reactions with shadow prices are:')
            for rid, sp in zip(self.kegg_model.rids, params['reaction prices'].flat):
                if sp:
                    logging.error('\t%s : %g' % (rid, sp))

            conc = list(map(np.exp, params['ln concentrations'].flat))
            for cid, sp, C in zip(self.kegg_model.cids, params['compound prices'].flat, conc):
                if sp and not self.cid2external[cid]:
                    logging.error('\t[%30s] : %5.1e < %5.1e < %5.1e M' %
                        (self.kegg2met[cid], self.cid2min_bound[cid], C, self.cid2max_bound[cid]))
            raise ThermodynamicallyInfeasibleError()
        return params['ln concentrations']

    def ECM(self, lnC0=None, n_iter=10):
        if lnC0 is None:
            lnC0 = self.MDF()
            logging.info('initializing ECM using MDF result')
        
        return self.ecf.ECM(lnC0, n_iter=n_iter)

    def ECF(self, lnC):
        return self.ecf.ECF(lnC)

    @staticmethod
    def _nanfloat(x):
        if type(x) == float:
            return x
        if type(x) == int:
            return float(x)
        if type(x) == str:
            if x.lower() in ['', 'nan']:
                return np.nan
            else:
                return float(x)
        else:
            raise ValueError('unrecognized type for value: ' + str(type(x)))

    @staticmethod
    def _MappingToCanonicalEnergyUnits(unit):
        """
            Assuming the canonical units for concentration are Molar

            Returns:
                A function that converts a single number or string to the
                canonical units
        """
        if unit == 'kJ/mol':
            return lambda x: ECMmodel._nanfloat(x)
        if unit == 'kcal/mol':
            return lambda x: ECMmodel._nanfloat(x)*4.184

        raise ValueError('Cannot convert these units to kJ/mol: ' + unit)

    @staticmethod
    def _MappingToCanonicalMolecularWeightUnits(unit):
        """
            Assuming the canonical units for concentration are Molar

            Returns:
                A function that converts a single number or string to the
                canonical units
        """
        if unit == 'Da':
            return lambda x: ECMmodel._nanfloat(x)
        if unit == 'kDa':
            return lambda x: ECMmodel._nanfloat(x)*1e3
        if unit == 'MDa':
            return lambda x: ECMmodel._nanfloat(x)*1e6
        if unit == 'mDa':
            return lambda x: ECMmodel._nanfloat(x)*1e-3

        raise ValueError('Cannot convert these units to M: ' + unit)

    @staticmethod
    def _MappingToCanonicalConcentrationUnits(unit):
        """
            Assuming the canonical units for concentration are Molar

            Returns:
                A function that converts a single number or string to the
                canonical units
        """
        if unit == 'M':
            return lambda x: ECMmodel._nanfloat(x)
        if unit == 'mM':
            return lambda x: ECMmodel._nanfloat(x)*1e-3
        if unit == 'uM':
            return lambda x: ECMmodel._nanfloat(x)*1e-6
        if unit == 'nM':
            return lambda x: ECMmodel._nanfloat(x)*1e-9

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
            return lambda x: ECMmodel._nanfloat(x)
        if unit == 'mM/s':
            return lambda x: ECMmodel._nanfloat(x)*1e-3
        if unit == 'mmol/gCDW/h':
            return lambda x: ECMmodel._nanfloat(x) / (CELL_VOL_PER_DW * 3600 * 1e3)
        if unit == 'mol/gCDW/h':
            return lambda x: ECMmodel._nanfloat(x) / (CELL_VOL_PER_DW * 3600)
        if unit == 'umol/gCDW/min':
            return lambda x: ECMmodel._nanfloat(x) / (CELL_VOL_PER_DW * 60 * 1e6)
        if unit == 'mmol/gCDW/min':
            return lambda x: ECMmodel._nanfloat(x) / (CELL_VOL_PER_DW * 60 * 1e3)

        raise ValueError('Cannot convert these units to M/s: ' + unit)

    def _GetMeasuredMetaboliteConcentrations(self):
        if self._validate_sbtab is None:
            raise Exception('cannot validate results because no validation data'
                            ' was given')
        unit = self._validate_sbtab.GetTableAttribute('Concentration', 'Unit')
        value_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(unit)

        # assume concentrations are in mM
        return self._validate_sbtab.GetDictFromTable(
            'Concentration', 'Compound:Identifiers:kegg.compound',
            'Value', value_mapping=value_mapping)

    def _GetMeasuredEnzymeConcentrations(self):
        if self._validate_sbtab is None:
            raise Exception('cannot validate results because no validation data'
                            ' was given')
        unit = self._validate_sbtab.GetTableAttribute('Concentration', 'Unit')
        value_mapping = ECMmodel._MappingToCanonicalConcentrationUnits(unit)

        return self._validate_sbtab.GetDictFromTable(
            'EnzymeConcentration', 'Reaction',
            'Value', value_mapping=value_mapping)

    def _GetVolumeDataForPlotting(self, lnC):
        enz_vols, met_vols = self.ecf.GetVolumes(lnC)

        enz_labels = list(map(self.kegg2rxn.get, self.kegg_model.rids))
        enz_data = sorted(zip(enz_vols.flat, enz_labels), reverse=True)
        enz_vols, enz_labels = zip(*enz_data)
        enz_colors = [(0.5, 0.8, 0.3)] * len(enz_vols)

        met_labels = list(map(self.kegg2met.get, self.kegg_model.cids))
        met_data = zip(met_vols.flat, met_labels)
        # remove H2O from the list and sort by descending volume
        met_data = sorted(filter(lambda x: x[1] != 'h2o', met_data))
        met_vols, met_labels = zip(*met_data)
        met_colors = [(0.3, 0.5, 0.8)] * len(met_vols)

        return (enz_vols + met_vols,
                enz_labels + met_labels,
                enz_colors + met_colors)

    def PlotVolumes(self, lnC, ax):
        width = 0.8
        vols, labels, colors = self._GetVolumeDataForPlotting(lnC)

        ax.bar(np.arange(len(vols)), vols, width, color=colors)
        ax.set_xticklabels(labels, size='medium', rotation=90)
        ax.set_ylabel('total weight [g/L]')

    def PlotVolumesPie(self, lnC, ax):
        vols, labels, colors = self._GetVolumeDataForPlotting(lnC)
        ax.pie(vols, labels=labels, colors=colors)
        ax.set_title('total weight [g/L]')

    def PlotEnzymeDemandBreakdown(self, lnC, ax, top_level=3, plot_measured=False):
        """
            A bar plot in log-scale showing the partitioning of cost between
            the levels of kinetic costs:
            1 - capacity
            2 - thermodynamics
            3 - saturation
            4 - allosteric
        """
        assert top_level in range(1, 5)

        costs = np.array(self.ecf.GetEnzymeCostPartitions(lnC))
        
        # give all reactions with zero cost a base value, which we will
        # also set as the bottom ylim, which will simulate a "minus infinity"
        # when we plot it in log-scale
        base = min(filter(None, costs[:, 0])) / 2.0
        idx_zero = (costs[:, 0] == 0)
        costs[idx_zero, 0] = base
        costs[idx_zero, 1:] = 1.0

        bottoms = np.hstack([np.ones((costs.shape[0], 1)) * base,
                             np.cumprod(costs, 1)])
        steps = np.diff(bottoms)

        labels = EnzymeCostFunction.ECF_LEVEL_NAMES[0:top_level]

        ind = range(costs.shape[0])    # the x locations for the groups
        width = 0.8
        ax.set_yscale('log')

        if plot_measured:
            all_labels = ['measured'] + labels
            meas_enz2conc = self._GetMeasuredEnzymeConcentrations()
            meas_conc = np.matrix(list(map(meas_enz2conc.get, self.kegg_model.rids)),
                                    dtype=float).T
            cmap = ColorMap(all_labels, saturation=0.7, value=1.0,
                            hues=[30.0/255, 170.0/255, 200.0/255, 5.0/255])
            ax.plot(ind, meas_conc, color=cmap['measured'], marker='o',
                    markersize=5, linewidth=0,
                    markeredgewidth=0.3, markeredgecolor=(0.3, 0.3, 0.3))
        else:
            all_labels = labels
            cmap = ColorMap(labels, saturation=0.7, value=0.8,
                            hues=[170.0/255, 200.0/255, 5.0/255])

        for i, label in enumerate(labels):
            ax.bar(ind, steps[:, i].flat, width,
                   bottom=bottoms[:, i].flat, color=cmap[label])

        ax.set_xticks(ind)
        xticks = list(map(self.kegg2rxn.get, self.kegg_model.rids))
        ax.set_xticklabels(xticks, size='medium', rotation=90)
        ax.legend(all_labels, loc='best', framealpha=0.2)
        ax.set_ylabel('enzyme demand [M]')
        ax.set_ylim(ymin=base)

    def ValidateMetaboliteConcentrations(self, lnC, ax, scale='log'):
        pred_conc = np.exp(lnC)

        meas_met2conc = self._GetMeasuredMetaboliteConcentrations()
        meas_conc = np.matrix(list(map(meas_met2conc.get, self.kegg_model.cids))).T

        # remove NaNs and zeros
        mask =  np.nan_to_num(meas_conc) > 0
        mask &= np.nan_to_num(pred_conc) > 0
        mask &= np.diff(self.ecf.lnC_bounds) > 1e-9 # remove compounds with fixed concentrations

        labels = list(map(self.kegg2met.get, self.kegg_model.cids))
        PlotCorrelation(ax, meas_conc, pred_conc, labels, mask, scale=scale)
        ax.set_xlabel('measured [M]')
        ax.set_ylabel('predicted [M]')

    def ValidateEnzymeConcentrations(self, lnC, ax, scale='log'):
        pred_conc = self.ecf.ECF(lnC)

        meas_enz2conc = self._GetMeasuredEnzymeConcentrations()
        meas_conc = np.matrix(list(map(meas_enz2conc.get, self.kegg_model.rids))).T

        labels = list(map(self.kegg2rxn.get, self.kegg_model.rids))
        PlotCorrelation(ax, meas_conc, pred_conc, labels, scale=scale)

        ax.set_xlabel('measured [M]')
        ax.set_ylabel('predicted [M]')

    def ToSBtab(self, lnC, filename, document_name=''):
        sys.path.append(os.path.expanduser('~/git/SBtab/python'))
        from SBtab import SBtabTable
        met_data = []
        for i, cid in enumerate(self.kegg_model.cids):
            met_name = self.kegg2met[cid]
            met_data.append(('concentration', met_name, cid, np.exp(lnC[i, 0])))
        met_df = pd.DataFrame(columns=['QuantityType', 'Compound',
                                        'Compound:Identifiers:kegg.compound', 'ecm'],
                               data=met_data)
        met_sbtab = SBtabTable.fromDataFrame(met_df,
            document_name,
            'Quantity',
            'Predicted concentrations',
            'ECM metabolic state',
            'M')
        met_sbtab.createDataset()
        
        enz_conc = self.ecf.ECF(lnC)
        enz_data = []
        for i, rid in enumerate(self.kegg_model.rids):
            rxn_name = self.kegg2rxn[rid]
            enz_data.append(('concentration of enzyme', rxn_name, rid, enz_conc[i, 0]))
        enz_df = pd.DataFrame(columns=['QuantityType', 'Reaction',
                                       'Reaction:Identifiers:kegg.reaction', 'ecm'],
                               data=enz_data)
        enz_sbtab = SBtabTable.fromDataFrame(enz_df,
            document_name,
            'Quantity',
            'Predicted enzyme levels',
            'ECM metabolic state',
            'M')
        enz_sbtab.createDataset()
        
        # TODO: write both results as a concatenated SBtab (i.e. one file 
        # with two tables)
        met_sbtab.writeSBtab('tsv', filename + 'met')
        enz_sbtab.writeSBtab('tsv', filename + 'enz')

    def WriteHtmlTables(self, lnC, html):
        meas_enz2conc = self._GetMeasuredEnzymeConcentrations()
        meas_conc = np.matrix(map(lambda r: meas_enz2conc.get(r, np.nan),
                               self.kegg_model.rids)).T
        data_mat = np.hstack([self.ecf.flux,
                              meas_conc,
                              self.ecf.ECF(lnC),
                              self.ecf._DrivingForce(lnC),
                              self.ecf.GetEnzymeCostPartitions(lnC)])

        data_mat[:, 0] *= 1e3 # convert flux from M/s to mM/s
        data_mat[:, 1] *= 1e6 # convert measured enzyme conc. from M to uM
        data_mat[:, 2] *= 1e6 # convert predicted enzyme conc. from M to uM
        data_mat[:, 4] *= 1e6 # convert capacity term from M to uM

        headers = ['reaction', 'KEGG ID', 'flux [mM/s]',
                   'measured enz. conc. [uM]',
                   'predicted enz. conc. [uM]',
                   'driving force [kJ/mol]', 'capacity [uM]'] + \
                   EnzymeCostFunction.ECF_LEVEL_NAMES[1:]
        values = zip(map(self.kegg2rxn.get, self.kegg_model.rids),
                     self.kegg_model.rids,
                     *data_mat.T.tolist())
        values.append(['total', '', '', data_mat[:, 1].sum(), data_mat[:, 2].sum(),
                       '', '', '', '', ''])

        rowdicst = [dict(zip(headers, v)) for v in values]
        html.write_table(rowdicst, headers=headers, decimal=3)

        meas_met2conc = self._GetMeasuredMetaboliteConcentrations()
        meas_conc = np.matrix(list(map(lambda r: meas_met2conc.get(r, np.nan),
                                    self.kegg_model.cids))).T
        data_mat = np.hstack([meas_conc,
                              np.exp(lnC),
                              np.exp(self.ecf.lnC_bounds)])

        data_mat *= 1e3 # convert all concentrations from M to mM

        headers = ['compound', 'KEGG ID', 'measured conc. [mM]',
                   'predicted conc. [mM]',
                   'lower bound [mM]', 'upper bound [mM]']
        values = zip(map(self.kegg2met.get, self.kegg_model.cids),
                     self.kegg_model.cids,
                     *data_mat.T.tolist())
        rowdicst = [dict(zip(headers, v)) for v in values]
        headers = ['compound', 'KEGG ID', 'measured conc. [mM]',
                   'lower bound [mM]', 'predicted conc. [mM]', 'upper bound [mM]']
        html.write_table(rowdicst, headers=headers)

