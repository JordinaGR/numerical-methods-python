import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def odef(t, y):
    return np.array((y[1], -y[0]))

def Euler(odef, initialTime, finalTime, alpha, nOfSteps):
    h = (finalTime-initialTime)/nOfSteps
    #t = np.arange(initialTime,finalTime+h,h)
    t = np.linspace(initialTime, finalTime, nOfSteps + 1)
    y = np.zeros((nOfSteps+1,len(alpha)))
    y[0,:] = alpha
    for i in np.arange(0,nOfSteps):
        y[i+1,:] = y[i,:] + h*odef(t[i], y[i,:]) 

    return t, y

def Heune(odef, initialTime, finalTime, alpha, nOfSteps):
    h = (finalTime-initialTime)/nOfSteps
    #t = np.arange(initialTime,finalTime+h,h)
    t = np.linspace(initialTime, finalTime, nOfSteps + 1)
    y = np.zeros((nOfSteps+1,len(alpha)))
    y[0,:] = alpha
    for i in np.arange(0,nOfSteps):
        y[i+1,:] = y[i,:] + (h/2)*(odef(t[i], y[i,:]) + odef(t[i+1], y[i, :] + h*odef(t[i], y[i,:]))) 

    return t, y

nOfSteps = 100
alpha = np.array((1., 0.))
initialTime = 0.
finalTime = 1.
t, y = Heune(odef, initialTime, finalTime, alpha, nOfSteps)

'''plt.plot(y[:, 0], y[:, 1], label = "euler")
plt.legend()
plt.show()
'''
print(y[-1, 0])  # ultima fila (ultima iteracio, interval) primera casella y:[0] y':[1]

# valor exacte
t_span = (0, 1)       # Des de t=0 fins a t=1
y0 = [1, 0]           # y(0)=1, y'(0)=0
solucio = solve_ivp(odef, t_span, y0, rtol=1e-13, atol=1e-13)
valor_final_scipy = solucio.y[0][-1]  # El valor de y a t=1

print(f"Solució amb SciPy (solve_ivp): {valor_final_scipy:.15f}")
