import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

'''
A = -1/3 * np.array([[4., 1.], [2., 5.]])
def odef(t,y):
    return A @ y
    #return np.array([y[1], -y[0]])
#    return np.array([y[0]+y[1], 2*y[0]*y[1]])
alpha = np.array([1,0]) #initial condition y(0)
finalTime=2
'''

A = np.array([[-4., -3.], [-6., -7.]])
def odef(t,y):
    return A @ y
    #return np.array([y[1], -y[0]])
#    return np.array([y[0]+y[1], 2*y[0]*y[1]])
alpha = np.array([1,1]) #initial condition y(0)
finalTime=2

#Euler method
nOfSteps=10
h=finalTime/nOfSteps
t=np.arange(0,finalTime+h,h)
y = np.zeros((nOfSteps+1,len(alpha)))
y[0,:] = alpha
for i in np.arange(0,nOfSteps):
    y[i+1,:] = y[i,:] + h*odef(t[i], y[i,:])    
# Plots
plt.plot(t,y[:,0],'*-')
plt.plot(t,y[:,1],'*-')
plt.legend(['y1','y2'])
plt.title('Euler method')
plt.show()
endError=y[-1,0]
print('Error final time = %0.1e' %(endError) )

# euler enrere cas lineal
y1 = np.zeros((nOfSteps+1, len(alpha)))
y1[0, :] = alpha
M = np.eye(len(alpha))-h * A
for i in np.arange(0, nOfSteps):
    y1[i+1, :] = np.linalg.solve(M, y1[i, :])


# #ODE solution with RK45
# sol = solve_ivp(odef, [0, finalTime], alpha,method='RK45') #,rtol=1.e-8)
# tref=sol.t
# yref=sol.y
# #Comparison plots
# plt.plot(t,y[:,0],'*-')
# plt.plot(tref,yref[0,:],'-o')
# plt.legend(['Euler','RK45'])
# plt.title('Comparison')
# plt.show()
