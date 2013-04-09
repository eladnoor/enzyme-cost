import cvxpy
import numpy as np

def main():
    N = 10
    S = cvxpy.matrix(np.zeros((N+1, N)))
    for i in xrange(N):
        S[i, i] = -1
        S[i+1, i] = 1
        
    d0 = cvxpy.matrix([-10, -10, -10, -10, 0, 0.5, 0.5, -10, -10, -10]).T
    
    lnC = cvxpy.variable(N+1, 1, name='lnC')
    
    d = d0 + S.T * lnC
    constraints = [cvxpy.geq(lnC, 0), cvxpy.leq(lnC, 1)]

    x = -cvxpy.log(1 - cvxpy.exp(-d))
    E = cvxpy.log_sum_exp(x)
    
    objective = cvxpy.minimize(E)
    
    program = cvxpy.program(objective, constraints)
    program.solve(quiet=True)
    
    print d.value
    print lnC.value
    
if __name__ == "__main__":
    main()