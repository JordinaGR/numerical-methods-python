# sistemes d’equacions no lineals, mètode de newton-raphson i newton-raphson modificat

import numpy as np
import matplotlib.pyplot as plt

def f(x):
    v = [
        6*x[0] - 2*np.cos(x[1]*x[2]) - 1,
        9*x[1] + np.sqrt(x[0]**2 + np.sin(x[2]) + 1.06) + 0.9,
        60*x[2] + 3*np.exp(x[0]*x[1]) + 10*np.pi - 3]
   
    return v

def jacobiana(x):
    J = np.zeros((3, 3))
   
    J[0][0] = 6
    J[0][1] = 2*np.sin(x[1]*x[2])*x[2]
    J[0][2] = 2*np.sin(x[1]*x[2])*x[1]
   
    J[1][0] = (1/np.sqrt(x[0]**2 + np.sin(x[2]) + 1.06))*x[0]
    J[1][1] = 9
    J[1][2] = (1/(2*np.sqrt(x[0]**2 + np.sin(x[2]) + 1.06)))*(np.cos(x[2]))
   
    J[2][0] = 3*np.exp(x[0]*x[1])*(x[1])
    J[2][1] = 3*np.exp(x[0]*x[1])*(x[0])
    J[2][2] = 60
   
    return J

def jacobianaNum(x, tol=1.8):
    J = np.zeros((3, 3))
    f0 = np.array(f(x))

    for j in range(3):
        x_mod = np.copy(x)
        x_mod[j] += tol
        f_mod = np.array(f(x_mod))
        J[:, j] = (f_mod - f0) / tol 

    return J

def NewtonRaphson(x0, iter, tol=1e-7):
    x1 = x0
    r = []
    
    for i in range(iter):
        if (np.linalg.norm(f(x0)) <= tol): break
        sol = np.linalg.solve(jacobiana(x0), f(x0))     
        x1 = x0 - sol
        r.append(np.linalg.norm(x1-x0) / np.linalg.norm(x1))
        x0 = x1
        
    return x1, r

def NewtonRaphsonNum(x0, iter, tol=1e-7):
    x1 = x0
    r = []
    
    for i in range(iter):
        if (np.linalg.norm(f(x0)) <= tol): break
        sol = np.linalg.solve(jacobianaNum(x0), f(x0))     
        x1 = x0 - sol
        r.append(np.linalg.norm(x1-x0) / np.linalg.norm(x1))
        x0 = x1
        
    return x1, r

def NewtonRaphsonMod(x0, iter, tol=1e-7):
    J = jacobiana(x0)
    x1 = x0
    r = []
    
    for i in range(iter):
        if (np.linalg.norm(f(x0)) <= tol): break
        sol = np.linalg.solve(J, f(x0))     
        x1 = x0 - sol
        r.append(np.linalg.norm(x1-x0) / np.linalg.norm(x1))
        x0 = x1
       
    return x1, r

x = [0,0,0]
niter = 50
sol, error = NewtonRaphson(x, 50)

# 3 les arrels als grafics sembla que estiguin al (0,0,0)
'''xpoints = np.linspace(-5, 5, 100)
arr = [xpoints, xpoints, xpoints]
plt.plot(xpoints, f(arr)[0])
plt.plot(xpoints, f(arr)[1])
plt.plot(xpoints, f(arr)[2])
plt.show()
'''

# 4 calcula l'error i grafic de l'error, l'ultim cas no convergeix
sol1, error1 = NewtonRaphson([0,0,0], niter)
sol2, error2 = NewtonRaphson([1, 1, 1], niter)
sol3, error3 = NewtonRaphson([5, 5, 5], niter)
#sol4, error4 = NewtonRaphson([-15, 15, -15], niter)
'''
# 5 newton raphson modificat
sol1m, error1m = NewtonRaphsonMod([0,0,0], niter)
sol2m, error2m = NewtonRaphsonMod([1, 1, 1], niter)
sol3m, error3m = NewtonRaphsonMod([5, 5, 5], niter)
#sol4m, error4m = NewtonRaphsonMod([-15, 15, -15], niter)
'''
# derivada numèrica
sol1m, error1m = NewtonRaphsonNum([0,0,0], niter)
sol2m, error2m = NewtonRaphsonNum([1, 1, 1], niter)
sol3m, error3m = NewtonRaphsonNum([5, 5, 5], niter)

plt.semilogy(range(1, len(error1) + 1), error1, marker='o')
plt.semilogy(range(1, len(error2) + 1), error2, marker='o')
plt.semilogy(range(1, len(error3) + 1), error3, marker='o')
#plt.semilogy(range(1, len(error4) + 1), error4, marker='o')
plt.show()

plt.semilogy(range(1, len(error1m) + 1), error1m, marker='o')
plt.semilogy(range(1, len(error2m) + 1), error2m, marker='o')
plt.semilogy(range(1, len(error3m) + 1), error3m, marker='o')
#plt.semilogy(range(1, len(error4) + 1), error4, marker='o')
plt.show()

# amb la derivada numèrica observem que calen més iteracions i la convergència sembla lineal. 
# amb la h inferior a 1 obtenim la matriu jacobiana singular i no es pot resoldre el sistema per fer servir el mètode de newton
