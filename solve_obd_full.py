import cvxpy
import numpy as np
import matplotlib.pyplot as plt
from optimized_bottleneck_driving_force import Pathway

def main():
    G0 = [-20, -20, 3, 6, -20, -20, 1, 1, 1, -20, -20, -20] # in units of RT
    N = len(G0)
    fluxes = [1] * N
    S = np.zeros((N+1, N))
    for i in xrange(N):
        S[i, i] = -1
        S[i+1, i] = 1
    x_min = np.ones((N+1, 1)) * -6 * np.log(10)
    x_max = np.ones((N+1, 1)) * -2 * np.log(10)
    path = Pathway(S, fluxes, G0, x_min, x_max)

    obd, params = path.FindOBD()

    cba, params = path.FindCBA()
    lnC_opt = params['concentrations']
    G_opt = list(params['Gibbs energies'].flat)

    print "OBD = %.3g" % obd
    print "OBD(min enz.) = %.3g" % -np.max(G_opt)
    return
    
    plt.subplot(121)
    plt.plot(np.cumsum([0] + G0), '-b', label='$\Delta_r G\'^\circ$')
    plt.plot(np.cumsum([0] + G_opt), '-g', label='$\Delta_r G\'$')
    plt.legend()
    plt.title('CBA = %2.g' % cba)
    plt.subplot(122)
    plt.plot(range(N+1), lnC_opt, '*r')
    plt.yscale('log')
    plt.ylabel('concentration [M]')
    plt.ylim(1e-7, 1e-1)
    plt.tight_layout(1.3)

    plt.show()
    
if __name__ == "__main__":
    main()