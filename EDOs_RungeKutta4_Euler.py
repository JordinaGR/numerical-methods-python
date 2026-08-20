import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

def odef(t, y):
    return np.array((y[1], 2*y[1] - 2*y[0]))

def f(x):
    return np.exp(x)*np.sin(x)

def Euler(odef, initialTime, finalTime, alpha, nOfSteps):
    h = (finalTime-initialTime)/nOfSteps
    t = np.arange(initialTime,finalTime+h,h)
    y = np.zeros((nOfSteps+1,len(alpha)))
    y[0,:] = alpha
    for i in np.arange(0,nOfSteps):
        y[i+1,:] = y[i,:] + h*odef(t[i], y[i,:]) 

    return t, y

def RK4(odef, x, y, h):
    k1 = odef(x, y)
    k2 = odef(x +h/2, y + (h/2)*k1)
    k3 = odef(x + h/2, y + (h/2)*k2)
    k4 = odef(x + h, y + h*k3)
    return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

def RK4_solve(odef, initialTime, finalTime, alpha, nOfSteps):
    h = (finalTime - initialTime) / nOfSteps
    t = np.linspace(initialTime, finalTime, nOfSteps+1)
    
    y = np.zeros((nOfSteps+1, len(alpha)))
    y[0,:] = alpha

    for i in range(nOfSteps):
        y[i+1,:] = RK4(odef, t[i], y[i,:], h)

    return t, y


def euler_enrere(odef, initialTime, finalTime, alpha, nOfSteps):
    h = (finalTime - initialTime) / nOfSteps
    t = np.linspace(initialTime, finalTime, nOfSteps+1)
    y = np.zeros((nOfSteps+1, len(alpha)))
    y[0,:] = alpha
    
    for i in range(nOfSteps):
        
        def F(z):
            return z - (y[i,:] + h * odef(t[i+1], z))
        
        z0 = y[i,:] + h * odef(t[i], y[i,:])
        y[i+1, :] = fsolve(F, z0)

    return t, y

nOfSteps = 10
alpha = np.array((0., 1.))
initialTime = 0.
finalTime = 1.
exctVal = f(finalTime)

# 10 passos
'''
t, y = Euler(odef, initialTime, finalTime, alpha, nOfSteps)
print(np.linalg.norm(exctVal-y[-1,:][0]))
xpoints = np.linspace(0, 1, 10)
plt.plot(xpoints, f(xpoints))
plt.plot(t, y[:, 0])
'''
# 20 passos
nOfSteps = 200
t, y = Euler(odef, initialTime, finalTime, alpha, nOfSteps)
print(np.linalg.norm(exctVal-y[-1,:][0]))
'''
xpoints = np.linspace(0, 1, 10)
plt.plot(xpoints, f(xpoints))
plt.plot(t, y[:, 0])
plt.show()
'''
errors = []
for i in range(1, 20, 1):
    t, y = Euler(odef, initialTime, finalTime, alpha, i)
    errors.append(np.linalg.norm(exctVal-y[-1, :][0]))
    
'''plt.loglog(range(1, len(errors) + 1), errors, marker='o')
plt.show()
'''
t, y = RK4_solve(odef, initialTime, finalTime, alpha, nOfSteps)
print(np.linalg.norm(exctVal-y[-1, :][0]))

t, y = euler_enrere(odef, initialTime, finalTime, alpha, nOfSteps)
print(np.linalg.norm(exctVal-y[-1, :][0]))
