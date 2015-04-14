import numpy as np
import pulp
from component_contribution.thermodynamic_constants import default_RT as RT

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

        flux_sign = map(np.sign, self.fluxes.flat)
        active_fluxes = np.abs(np.matrix(flux_sign)).T
        I_dir = np.matrix(np.diag(flux_sign))
        
        A = np.matrix(np.vstack([np.hstack([I_dir * self.S.T, active_fluxes]),
                                 np.hstack([np.eye(self.Nc),  np.zeros((self.Nc, 1))]),
                                 np.hstack([-np.eye(self.Nc), np.zeros((self.Nc, 1))])]))
        b = np.matrix(np.vstack([-(I_dir * self.dG0) / RT,
                                  self.lnC_bounds[:, 1:],
                                 -self.lnC_bounds[:, :1]]))
        c = np.matrix(np.vstack([np.zeros((self.Nc, 1)),
                                 np.ones((1, 1))]))
       
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
        _l = pulp.LpVariable.dicts("lnC", ["%d" % i for i in xrange(self.Nc)])
        B = pulp.LpVariable("B")
        lnC = [_l["%d" % i] for i in xrange(self.Nc)] + [B]
        
        for j in xrange(A.shape[0]):
            row = [A[j, i] * lnC[i] for i in xrange(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        objective = pulp.lpSum([c[i] * lnC[i] for i in xrange(A.shape[1])])
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
                                  ["%d" % i for i in xrange(self.Nr)],
                                  lowBound=0)

        z = pulp.LpVariable.dicts("z", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)

        u = pulp.LpVariable.dicts("u", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)
        
        y = [w["%d" % i] for i in xrange(self.Nr)] + \
            [z["%d" % i] for i in xrange(self.Nc)] + \
            [u["%d" % i] for i in xrange(self.Nc)]
        
        for i in xrange(A.shape[1]):
            row = [A[j, i] * y[j] for j in xrange(A.shape[0])]
            lp += (pulp.lpSum(row) == c[i, 0]), "dual_%02d" % i

        objective = pulp.lpSum([b[i] * y[i] for i in xrange(A.shape[0])])
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
        lp_primal.solve(pulp.CPLEX(msg=0))
        if lp_primal.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF primal")
            
        mdf = pulp.value(lnC[-1])
        lnC = np.matrix(map(pulp.value, lnC[0:self.Nc])).T
    
        lp_dual, w, z, u = self._MakeMDFProblemDual()
        lp_dual.solve(pulp.CPLEX(msg=0))
        if lp_dual.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF dual")
        reaction_prices = np.matrix([pulp.value(w["%d" % i]) for i in xrange(self.Nr)]).T
        compound_prices = np.matrix([pulp.value(z["%d" % j]) for j in xrange(self.Nc)]).T - \
                          np.matrix([pulp.value(u["%d" % j]) for j in xrange(self.Nc)]).T
        
        params = {'MDF': mdf,
                  'ln concentrations' : lnC,
                  'reaction prices' : reaction_prices,
                  'compound prices' : compound_prices}
        return mdf, params
