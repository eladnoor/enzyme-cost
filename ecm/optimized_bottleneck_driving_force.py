import numpy as np
import pulp
from ecm.util import RT

class Pathway(object):

    def __init__(self, S, fluxes, dG0, lnC_bounds):
        """
            All inputs should be of type numpy.array:

            S            - the stoichiometric matrix
            fluxes       - the relative flux in each reaction
            G0           - the standard reaction Gibbs energies (in units of RT)
            lnC_bounds   - the natural log of the lower and upper bounds on
                           compounds concentrations (relative to 1M)
        """
        self.Nc, self.Nr = S.shape
        assert fluxes.shape     == (self.Nr, 1)
        assert dG0.shape        == (self.Nr, 1)
        assert lnC_bounds.shape == (self.Nc, 2)

        self.S = S
        self.fluxes = fluxes
        self.dG0 = dG0
        self.lnC_bounds = lnC_bounds

    def _MakeDrivingForceConstraints(self):
        """
            Generates the A matrix and b & c vectors that can be used in a
            standard form linear problem:
                max          c'x
                subject to   Ax >= b
                             x >= 0
        """
        # we need to take special care for reactions whose flux is 0.
        # since we don't want them to constrain the MDF to be 0.

        flux_sign = list(map(np.sign, self.fluxes.flat))
        active_fluxes = np.abs(np.array(flux_sign, ndmin=2)).T
        I_dir = np.diag(flux_sign)

        A = np.vstack([np.hstack([np.dot(I_dir, self.S.T), active_fluxes]),
                       np.hstack([np.eye(self.Nc),  np.zeros((self.Nc, 1))]),
                       np.hstack([-np.eye(self.Nc), np.zeros((self.Nc, 1))])])
        b = np.vstack([-(np.dot(I_dir, self.dG0)) / RT,
                       self.lnC_bounds[:, 1:],
                       -self.lnC_bounds[:, :1]])
        c = np.vstack([np.zeros((self.Nc, 1)), np.ones((1, 1))])

        return A, b, c

    def _MakeMDFProblem(self):
        """Create a PuLP problem for finding the Maximal Thermodynamic
        Driving Force (MDF).

        Does not set the objective function... leaves that to the caller.

        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.

        Returns:
            A tuple (dgf_var, motive_force_var, problem_object).
        """
        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints()

        lp = pulp.LpProblem("MDF", pulp.LpMaximize)

        # ln-concentration variables
        _l = pulp.LpVariable.dicts("lnC", ["%d" % i for i in range(self.Nc)])
        B = pulp.LpVariable("B")
        lnC = [_l["%d" % i] for i in range(self.Nc)] + [B]

        for j in range(A.shape[0]):
            row = [A[j, i] * lnC[i] for i in range(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j

        objective = pulp.lpSum([c[i] * lnC[i] for i in range(A.shape[1])])
        lp.setObjective(objective)

        #lp.writeLP("res/MDF_primal.lp")

        return lp, lnC

    def _MakeMDFProblemDual(self):
        """Create a CVXOPT problem for finding the Maximal Thermodynamic
        Driving Force (MDF).

        Does not set the objective function... leaves that to the caller.

        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.

        Returns:
            A tuple (dgf_var, motive_force_var, problem_object).
        """
        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints()

        lp = pulp.LpProblem("MDF", pulp.LpMinimize)

        w = pulp.LpVariable.dicts("w",
                                  ["%d" % i for i in range(self.Nr)],
                                  lowBound=0)

        z = pulp.LpVariable.dicts("z",
                                  ["%d" % i for i in range(self.Nc)],
                                  lowBound=0)

        u = pulp.LpVariable.dicts("u",
                                  ["%d" % i for i in range(self.Nc)],
                                  lowBound=0)

        y = [w["%d" % i] for i in range(self.Nr)] + \
            [z["%d" % i] for i in range(self.Nc)] + \
            [u["%d" % i] for i in range(self.Nc)]

        for i in range(A.shape[1]):
            row = [A[j, i] * y[j] for j in range(A.shape[0])]
            lp += (pulp.lpSum(row) == c[i, 0]), "dual_%02d" % i

        objective = pulp.lpSum([b[i] * y[i] for i in range(A.shape[0])])
        lp.setObjective(objective)

        #lp.writeLP("res/MDF_dual.lp")

        return lp, w, z, u

    def FindMDF(self):
        """
            Find the MDF (Optimized Bottleneck Driving-Force)

            Returns:
                A pair of the resulting MDF (in units of RT)
                and a dictionary of all the
                parameters and the resulting MDF value
        """
        lp_primal, lnC = self._MakeMDFProblem()
        lp_primal.solve()
        lp_primal.writeLP('/tmp/mdf.lp')
        if lp_primal.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF primal")

        mdf = pulp.value(lnC[-1])
        lnC = np.array(list(map(pulp.value, lnC[0:self.Nc])), ndmin=2).T

        lp_dual, w, z, u = self._MakeMDFProblemDual()
        lp_dual.solve()
        if lp_dual.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF dual")
        reaction_prices = np.array([pulp.value(w["%d" % i]) for i in range(self.Nr)], ndmin=2).T
        compound_prices = np.array([pulp.value(z["%d" % j]) for j in range(self.Nc)], ndmin=2).T - \
                          np.array([pulp.value(u["%d" % j]) for j in range(self.Nc)], ndmin=2).T

        params = {'MDF': mdf,
                  'ln concentrations' : lnC,
                  'reaction prices' : reaction_prices,
                  'compound prices' : compound_prices}
        return mdf, params
