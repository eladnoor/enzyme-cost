import numpy as np
import matplotlib.pyplot as plt
from ecm.optimized_bottleneck_driving_force import Pathway

def compare():
    G0 = np.random.randn(10, 1)*5 - 3
#    G0 = [-20, -20, 3, 6, -20, -20, 1, 1, 1, -20, -20, -20] # in units of RT
    N = G0.shape[0]
    fluxes = [1] * N
    S = np.zeros((N+1, N))
    for i in range(N):
        S[i, i] = -1
        S[i+1, i] = 1
    x_min = np.ones((N+1, 1)) * -6 * np.log(10)
    x_max = np.ones((N+1, 1)) * -2 * np.log(10)
    lnC_bounds = np.hstack([x_min, x_max])
    path = Pathway(S, fluxes, G0, lnC_bounds)

    obd, params = path.FindOBD()
    if obd < 0:
        return None, None

    cba, params = path.FindCBA()
    lnC_opt = params['concentrations']
    G_opt = list(params['Gibbs energies'].flat)

    obd_cost = -np.max(G_opt)
    print(obd, obd_cost)
    return obd, obd_cost
    
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
    
def multi_compare():
    data = []
    for i in range(20):
        obd, obd_cost = compare()
        if obd is not None:
            data.append((obd, obd_cost))
    data = np.array(data)
    plt.plot(data[:,0], data[:,1], '.')
    plt.xlabel('OBD')
    plt.ylabel('bottleneck DF after enzyme cost optimization')
    plt.show()
    
def try_ECF():
#    G0 = np.random.randn(10, 1)*5 - 3
    G0 = np.array([-20, -20, -3, -6, -20, -20, -1, -1, -1, -20, -20, -20]).T # in units of RT
    N = G0.shape[0]
    fluxes = np.array([1] * N, ndmin=2).T
    S = np.zeros((N+1, N))
    for i in range(N):
        S[i, i] = -1
        S[i+1, i] = 1
    x_min = np.ones((N+1, 1)) * -6 * np.log(10)
    x_max = np.ones((N+1, 1)) * -2 * np.log(10)
    lnC_bounds = np.hstack([x_min, x_max])
    path = Pathway(S, fluxes, G0, lnC_bounds)
    mdf, params = path.FindMDF()
    print(mdf)
    
if __name__ == "__main__":
    #multi_compare()
    try_ECF()
