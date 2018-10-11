#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 11:06:07 2017

@author: noore
"""

import re
import numpy as np
from equilibrator_api import Reaction, ComponentContribution

class KeggParseException(Exception):
    pass

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
    for cid, count in ParseReactionFormulaSide(left).items():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count

    for cid, count in ParseReactionFormulaSide(right).items():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count

    return sparse_reaction

class KeggReaction(object):

    def __init__(self, sparse, arrow='<=>', rid=None):
        for cid, coeff in sparse.items():
            if not (isinstance(coeff, float) or isinstance(coeff, int)):
                raise ValueError('All values in KeggReaction must be integers or floats')
        self.sparse = {k: v for (k,v) in sparse.items() if v != 0}
        self.arrow = arrow
        self.rid = rid

    def keys(self):
        return self.sparse.keys()

    def dense(self, cids):
        s = np.zeros((len(cids), 1))
        for cid, coeff in self.items():
            s[cids.index(cid), 0] = coeff
        return s

    def items(self):
        return self.sparse.items()

class KeggModel(object):

    def __init__(self, S, cids, rids=None):
        self.S = S
        self.cids = cids
        self.rids = rids
        assert len(self.cids) == self.S.shape[0]
        if self.rids is not None:
            assert len(self.rids) == self.S.shape[1]

        # remove H+ from the stoichiometric matrix if it exists
        if 'C00080' in self.cids:
            i = self.cids.index('C00080')
            self.S = np.vstack((self.S[:i,:], self.S[i+1:,:]))
            self.cids.pop(i)

        self.reactions = []
        for i in range(self.S.shape[1]):
            kegg_id_to_coeff = {self.cids[j]: self.S[j, i]
                                for j in range(self.S.shape[0])
                                if self.S[j, i] != 0}
            self.reactions.append(Reaction(kegg_id_to_coeff))

        # remove H2O from the stoichiometric matrix if it exists
        if 'C00001' in self.cids:
            i = self.cids.index('C00001')
            self.S = np.vstack((self.S[:i,:], self.S[i+1:,:]))
            self.cids.pop(i)

    @staticmethod
    def from_kegg_reactions(kegg_reactions, has_reaction_ids=False):
        if has_reaction_ids:
            rids = [r.rid for r in kegg_reactions]
        else:
            rids = None

        cids = set()
        for reaction in kegg_reactions:
            cids = cids.union(reaction.keys())

        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        cids = sorted(cids)
        S = np.zeros((len(cids), len(kegg_reactions)))
        for i, reaction in enumerate(kegg_reactions):
            S[:, i:i+1] = reaction.dense(cids)

        return KeggModel(S, cids, rids)

    def get_rid2dG0(self, pH=7.5, ionic_strength=0.1):
        equilibrator = ComponentContribution(pH=pH, ionic_strength=ionic_strength)
        rid2dG0 = dict(zip(self.rids,
                           map(lambda r: equilibrator.dG0_prime(r)[0], self.reactions)))
        return rid2dG0
