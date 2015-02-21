import numpy as np
import pulp
import types
import scipy.optimize

def CastToColumnVector(v):
    if type(v) == types.ListType:
        return np.array(v, ndmin=2).T
    if type(v) in [np.ndarray, np.matrix]:
        return np.array(v.flat, ndmin=2).T
    else:
        raise ValueError('Can only cast lists or numpy arrays, not ' + str(type(v)))
    
class Pathway(object):
    
    def __init__(self, S, fluxes, G0, x_min, x_max):
        """
            All inputs should be of type numpy.array:
            
            S            - the stoichiometric matrix
            fluxes       - the relative flux in each reaction
            G0           - the standard reaction Gibbs energies (in units of RT)
            x_min, x_max - the natural log of the lower and upper bounds on 
                           compounds concentrations (relative to 1M)
        """
        self.S = S
        self.Nc = self.S.shape[0]
        self.Nr = self.S.shape[1]

        self.fluxes = CastToColumnVector(fluxes)
        self.G0 = CastToColumnVector(G0)
        self.x_min = CastToColumnVector(x_min)
        self.x_max = CastToColumnVector(x_max)

        assert self.fluxes.shape == (self.Nr, 1)
        assert self.G0.shape == (self.Nr, 1)
        assert self.x_min.shape == (self.Nc, 1)
        assert self.x_max.shape == (self.Nc, 1)

    def _MakeDrivingForceConstraints(self):
        """
            Generates the A matrix and b & c vectors that can be used in a 
            standard form linear problem:
                max          c'x
                subject to   Ax >= b
                             x >= 0
        """
        I_dir = np.matrix(np.diag([np.sign(x) for x in self.fluxes.flat]))
       
        A = np.matrix(np.vstack([np.hstack([I_dir * self.S.T, np.ones((self.Nr, 1)) ]),
                                 np.hstack([np.eye(self.Nc),  np.zeros((self.Nc, 1))]),
                                 np.hstack([-np.eye(self.Nc), np.zeros((self.Nc, 1))])]))
        b = np.matrix(np.vstack([-I_dir * self.G0,
                                 self.x_max,
                                 -self.x_min]))
        c = np.matrix(np.vstack([np.zeros((self.Nc, 1)),
                                 np.ones((1, 1))]))
       
        return A, b, c

    def _MakeOBEProblem(self):
        """Create a PuLP problem for finding the Maximal Thermodynamic
        Driving Force (OBE).
       
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
        l = pulp.LpVariable.dicts("l", ["%d" % i for i in xrange(self.Nc)])
        B = pulp.LpVariable("B")
        x = [l["%d" % i] for i in xrange(self.Nc)] + [B]
        
        for j in xrange(A.shape[0]):
            row = [A[j, i] * x[i] for i in xrange(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        objective = pulp.lpSum([c[i] * x[i] for i in xrange(A.shape[1])])
        lp.setObjective(objective)
        
        #lp.writeLP("../res/obe_primal.lp")
        
        return lp, x

    def _MakeOBEProblemDual(self):
        """Create a CVXOPT problem for finding the Maximal Thermodynamic
        Driving Force (OBE).
       
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
       
        lp = pulp.LpProblem("OBE", pulp.LpMinimize)
        
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
        
        #lp.writeLP("../res/obe_dual.lp")
        
        return lp, w, z, u

    def _GetTotalEnergyProblem(self, min_driving_force=0, objective=pulp.LpMinimize):
        # Create the driving force variable and add the relevant constraints
        A, b, _c = self._MakeDrivingForceConstraints()
       
        lp = pulp.LpProblem("OBE", objective)
        
        # ln-concentration variables
        l = pulp.LpVariable.dicts("l", ["%d" % i for i in xrange(self.Nc)])
        x = [l["%d" % i] for i in xrange(self.Nc)] + [min_driving_force]
        
        total_g = pulp.LpVariable("g_tot")
        
        for j in xrange(A.shape[0]):
            row = [A[j, i] * x[i] for i in xrange(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        total_g0 = float(np.dot(self.G0.T, self.fluxes))
        total_reaction = np.dot(self.S, self.fluxes)
        row = [total_reaction[i, 0] * x[i] for i in xrange(self.Nc)]
        lp += (total_g == total_g0 + pulp.lpSum(row)), "Total G"

        lp.setObjective(total_g)
        
        #lp.writeLP("../res/total_g.lp")
        
        return lp, total_g
        
    def FindMDF(self):
        """
            Find the MDF (Optimized Bottleneck Driving-Force)
       
            Returns:
                A pair of the resulting MDF (in units of RT)
                and a dictionary of all the
                parameters and the resulting MDF value
        """
        lp_primal, x = self._MakeOBEProblem()
        lp_primal.writeLP('res/mdf.lp')
        lp_primal.solve(pulp.CPLEX(msg=0))
        if lp_primal.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF primal")
            
        mdf = pulp.value(x[-1])
        conc = np.matrix([np.exp(pulp.value(x[j])) for j in xrange(self.Nc)]).T
    
        lp_dual, w, z, u = self._MakeOBEProblemDual()
        lp_dual.solve(pulp.CPLEX(msg=0))
        if lp_dual.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve OBE dual")
        reaction_prices = np.matrix([pulp.value(w["%d" % i]) for i in xrange(self.Nr)]).T
        compound_prices = np.matrix([pulp.value(z["%d" % j]) for j in xrange(self.Nc)]).T - \
                          np.matrix([pulp.value(u["%d" % j]) for j in xrange(self.Nc)]).T
        
        # find the maximum and minimum total Gibbs energy of the pathway,
        # under the constraint that the driving force of each reaction is >= OBE
        lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-6, pulp.LpMinimize)
        lp_total.solve(pulp.CPLEX(msg=0))
        if lp_total.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve total delta-G problem")
        min_tot_dg = pulp.value(total_dg)
    
        lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-6, pulp.LpMaximize)
        lp_total.solve(pulp.CPLEX(msg=0))
        if lp_total.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve total delta-G problem")
        max_tot_dg = pulp.value(total_dg)
        
        params = {'MDF': mdf,
                  'concentrations' : conc,
                  'reaction prices' : reaction_prices,
                  'compound prices' : compound_prices,
                  'maximum total dG' : max_tot_dg,
                  'minimum total dG' : min_tot_dg}
        return mdf, params
    
    def FindCBA(self):
        import cvxpy
        """
            Solves the minimal enzyme cost function assuming the rate law: 
            
            v_i = E_i * (1 - exp(G_i))
            
            where G is the reaction Gibbs energy (in units of RT) 
        """
        lnC = cvxpy.variable(self.Nc, 1, name='lnC')
        G = cvxpy.matrix(self.G0) + cvxpy.matrix(self.S).T * lnC
        x = cvxpy.exp(-cvxpy.log(1 - cvxpy.exp(G)))
        E = x.T * cvxpy.matrix(self.fluxes)
    
        constraints = [cvxpy.geq(lnC, cvxpy.matrix(self.x_min))] + \
                      [cvxpy.leq(lnC, cvxpy.matrix(self.x_max))] + \
                      [cvxpy.leq(G, 0.0)]
        
        objective = cvxpy.minimize(E)
        
        program = cvxpy.program(objective, constraints)
        program.solve(quiet=True)
        
        G_opt = list(G.value.flat)
        lnC_opt = list(lnC.value.flat)

        cba = E.value
        params = {'CBA': cba,
                  'concentrations' : np.exp(lnC.value),
                  'Gibbs energies' : G.value}
        return cba, params
    
    def calc_E(self, x):
        S = np.matrix(self.S)
        x = np.matrix(x.flat)
        thermo = 1 - np.exp(self.G0 + S.T * x.T)
        E = np.dot((1.0 / thermo).T, self.fluxes)
        return float(E)
    
    def FindECF(self):
        """
            Solves the minimal enzyme cost function assuming the rate law: 
            
            v_i = E_i * (1 - exp(G_i)) * Pkin(x)
            
            where 'G' is the reaction Gibbs energy (in units of RT) 
            and 'Pkin' is a polynomial representing the kinetic term 
            where 'x' are the log-concentrations of metabolites
        """
        
        fun = lambda x : self.calc_E(x)
        x0 = (self.x_min + self.x_max) / 2

        print fun(x0)

        res = scipy.optimize.minimize(fun, x0)
        if res.success:
            return res.x
        else:
            print "No solution found: " + res.message
            return None
