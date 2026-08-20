import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def initialCondition(x):
    return 2*np.sin(np.pi*x)

def buildSystem(a, b, nOfIntervals, nu, BCtype, BCvalue_a, BCvalue_b):
    Dx = (b - a) / nOfIntervals
    if BCtype == 'DD':
        dm = nOfIntervals - 1
    elif BCtype == 'NN':
        dm = nOfIntervals + 1
    A = (nu/Dx**2) * diags([1, -2, 1], [-1, 0, 1], shape=(dm,dm)).toarray()
    F = np.zeros(dm)
    if BCtype == 'DD':
        F[0] = nu * BCvalue_a / Dx**2
        F[-1] = nu * BCvalue_b / Dx ** 2
    elif BCtype == 'NN':
        A[0, 1] = 2 * nu / Dx**2
        A[-1, -2] = 2 * nu / Dx**2
        F[0] = 2 * nu * BCvalue_a / Dx
        F[-1] = 2 * nu * BCvalue_b / Dx
    return A, F


def fwEuler(A, F, U0, T, nOfTimeSteps, BCtype, BCvalue_a, BCvalue_b):
    dm = len(U0)
    U = np.zeros((dm, nOfTimeSteps+1))
    U[:, 0] = U0
    if BCtype == 'DD':
        U[0, :] = BCvalue_a
        U[-1, :] = BCvalue_b
        ind = np.arange(1, dm-1)
    elif BCtype == 'NN':
        ind = np.arange(dm)

    Dt = T / nOfTimeSteps
    for n in range(nOfTimeSteps):
        Un = U[ind, n]
        U[ind, n+1] = Un + Dt*(A@Un + F)
    return U

def bwEuler(A, F, U0, T, nOfTimeSteps, BCtype, BCvalue_a, BCvalue_b):
    dm = len(U0)
    U = np.zeros((dm, nOfTimeSteps+1))
    U[:, 0] = U0
    if BCtype == 'DD':
        U[0, :] = BCvalue_a
        U[-1, :] = BCvalue_b
        ind = np.arange(1, dm-1)
    elif BCtype == 'NN':
        ind = np.arange(dm)

    Dt = T / nOfTimeSteps
    for n in range(nOfTimeSteps):
        Un = U[ind, n]
        K = np.eye(len(A[0, :])) -Dt*A
        b = Un + Dt*F
        U[ind, n+1] = np.linalg.solve(K,b)
    return U

# Definicio del problema
a = 0
b = 1
finalT = 0.4
nu = 1
BC = 'DD'
ua = 0
ub = 0
# Discretitzacio
nOfIntervals = 10
nOfTimeSteps = 40

A, F = buildSystem(a, b, nOfIntervals, nu, 'DD', ua, ub)
x = np.linspace(a, b, nOfIntervals+1)
U0 = initialCondition(x)
U_fw = fwEuler(A, F, U0, finalT, nOfTimeSteps, 'DD', ua, ub)
U_bw = bwEuler(A, F, U0, finalT, nOfTimeSteps, 'DD', ua, ub)


# observem que amb euler andavant no convergeix pero amb euler enrrere si
'''
plt.plot(x, U_fw)
plt.show()

plt.plot(x, U_bw)
plt.show()
'''
# per euler endavant veiem que At/Ax^2 = 1 <! 1/2 per tant no és estable
# en canvi per euler enrrere el mètode es incondicionalment estable

# 4
finalT = 0.1
nu = 0.25
U_bw = bwEuler(A, F, U0, finalT, nOfTimeSteps, 'DD', ua, ub)
plt.plot(x, U_bw)
plt.show()

finalT = 0.1
nu = 1
U_bw = bwEuler(A, F, U0, finalT, nOfTimeSteps, 'DD', ua, ub)
plt.plot(x, U_bw)
plt.show()

finalT = 0.1
nu = 4
U_bw = bwEuler(A, F, U0, finalT, nOfTimeSteps, 'DD', ua, ub)
plt.plot(x, U_bw)
plt.show()


nu = 1
finalT = 0.1
nOfIntervals = 10
mval = np.array([10*(2**m) for m in range (7)])
for m in mval:
    A, F = buildSystem(a, b, m, nu, 'DD', ua, ub)
    x = np.linspace(a, b, m+1)
    U0 = initialCondition(x)
    U_fw = fwEuler(A, F, U0, finalT, m, 'DD', ua, ub)
    U_bw = bwEuler(A, F, U0, finalT, m, 'DD', ua, ub)

    print(U_bw[5][-1])

    plt.plot(x, U_fw)
    plt.grid(True)
    plt.title('Euler')
    plt.show()
    
    plt.plot(x, U_bw)
    plt.grid(True)
    plt.title('Euler enrere')
    plt.show()
