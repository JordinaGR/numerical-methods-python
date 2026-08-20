# resolució d'equacions no lineals amb els mètodes de Newton, whittaker i mètode de la secant

import numpy as np
import matplotlib.pyplot as plt

# funcions de R en R

def f(x):
    return x**5 - 2*x**4 - 6*x**3 + 12*x**2 + 9*x - 18

def derivadaf(x):
    return 5*x**4 - 8*x**3 - 18*x**2 + 24*x + 9

def f2(x):
    return x**5 - 4*x**4 + 7*x**3 - 21*x**2 + 6*x + 18

def derivadaf2(x):
    return 5*x**4 - 16*x**3 + 21*x**2 - 42*x + 6

def errorRelatiu(a, b):
    return np.abs((a-b)/b)

def newton_iter(valorini, niter):
    error = []
    xanterior = valorini
   
    for i in range(niter+1):
        xseguent = xanterior - f(xanterior)/derivadaf(xanterior)      
        error.append(errorRelatiu(xanterior, xseguent))
        xanterior = xseguent
       
    return xseguent, error

def newton(valorini, niter, tol_x=1e-4, tol_f=1e-4): # metode de newton amb break tolerancia
    residus = []
    error = []
    xanterior = valorini
   
    for i in range(niter+1):
        if (abs(f2(xanterior)) <= tol_f):
            break
        residus.append(abs(f2(xanterior)))
        xseguent = xanterior - f2(xanterior)/derivadaf2(xanterior)      

        r_k = errorRelatiu(xanterior, xseguent)

        if (abs(r_k) <= tol_x):
            break

        error.append(r_k)
        xanterior = xseguent
       
    return xseguent, error, residus

def whittaker(valorini, niter, m):
    xanterior = valorini
    error = []
    for i in range(niter+1):
        xseguent = xanterior - f2(xanterior) / m
        error.append(errorRelatiu(xanterior, xseguent))
        xanterior = xseguent
       
    return xseguent, error

def secant(valorini, niter): # canviar segons convingui la f per f2
   
    error = []
    xanterior = valorini
    x = 0
   
    for i in range(niter+1):
        if (f2(xanterior)-f2(x) != 0): xseguent = xanterior - f2(xanterior)*((xanterior-x)/(f2(xanterior)-f2(x)))
        else: break
        error.append(errorRelatiu(xanterior, xseguent))
        x = xanterior
        xanterior = xseguent
        
    return xseguent, error

def newton_derivada(valorini, niter, h): # canviar segons convingui la f per f2
    error = []
    xanterior = valorini
   
    for i in range(niter+1):
        xseguent = xanterior - f2(xanterior)/((f2(xanterior+h)-f2(xanterior))/h)
        error.append(errorRelatiu(xanterior, xseguent))
        xanterior = xseguent
       
    return xseguent, error

niter = 20

# 2 dibuixa la gràfica de la funció, sembla que -1.5, 1, 3 és una aproximació inicial raonable
#xpoints = np.linspace(-1, 4, 100)
#plt.plot(xpoints, f2(xpoints))
#plt.show()

# 3 i 4 fes servir mètode de newton per trobar les arrels de (2) amb les seguents aproximacions
'''res1, e1 = newton_iter(-1, niter) # sol = -0,69
print(res1, f2(res1))

res2, e2 = newton_iter(2, niter) # sol = 1,23
print(res2, f2(res2))

res3, e3 = newton_iter(3, niter) # sol = 3,46
print(res3, f2(res3))

res4, e4 = newton_iter(2.5, niter) # sol = -0,69
print(res4, f2(res4))

plt.semilogy(range(1, len(e1) + 1), e1, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e2) + 1), e2, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e3) + 1), e3, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e4) + 1), e4, marker='o') # convergència quadràtica
plt.show()'''

# 6 
'''res1, e1, residus1 = newton(-1, niter) # sol = -0,69
print(res1, f2(res1))

res2, e2, residus2 = newton(2, niter) # sol = 1,23
print(res2, f2(res2))

res3, e3, residus3 = newton(3, niter) # sol = 3,46
print(res3, f2(res3))

res4, e4, residus4 = newton(2.5, niter) # sol = -0,69
print(res4, f2(res4))

plt.semilogy(range(1, len(e1) + 1), e1, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e2) + 1), e2, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e3) + 1), e3, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e4) + 1), e4, marker='o') # convergència quadràtica
plt.show()

plt.semilogy(range(1, len(residus1) + 1), residus1, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(residus2) + 1), residus2, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(residus3) + 1), residus3, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(residus4) + 1), residus4, marker='o') # convergència quadràtica
plt.show()'''

# ex addicionals 1
'''
res1, e1 = newton_iter(-1, niter) # sol = -1,73
print(res1, f(res1))

res2, e2 = newton_iter(1, niter) # sol = 1.73
print(res2, f(res2))

res3, e3 = newton_iter(3, niter) # sol = 2
print(res3, f(res3))

plt.semilogy(range(1, len(e1) + 1), e1, marker='o') # convergència lineal per tant arrel múltiple
plt.semilogy(range(1, len(e2) + 1), e2, marker='o') # convergència lineal per tant arrel múltiple
plt.semilogy(range(1, len(e3) + 1), e3, marker='o') # convergència quadràtica
plt.show()
'''

# addicionals 2 whittaker el millor valor de m es el que més s'asembla a la derivada
'''res1, e1 = whittaker(1, niter, -26)
print(res1, f2(res1))

res2, e2 = whittaker(1, niter, -32)
print(res2, f2(res2))

res3, e3 = whittaker(1, niter, -20)
print(res3, f2(res3))

res4, e4 = whittaker(1, niter, 150)
print(res4, f2(res4))

plt.semilogy(range(1, len(e1) + 1), e1, marker='o')
plt.semilogy(range(1, len(e2) + 1), e2, marker='o')
plt.semilogy(range(1, len(e3) + 1), e3, marker='o')
plt.semilogy(range(1, len(e4) + 1), e4, marker='o')
plt.show()'''

'''
res1, e1 = whittaker(2, niter, -26)
print(res1, f2(res1))

res2, e2 = whittaker(2, niter, -32)
print(res2, f2(res2))

res3, e3 = whittaker(2, niter, -20)
print(res3, f2(res3))

res4, e4 = whittaker(2, niter, 150)
print(res4, f2(res4))

plt.semilogy(range(1, len(e1) + 1), e1, marker='o')
plt.semilogy(range(1, len(e2) + 1), e2, marker='o')
plt.semilogy(range(1, len(e3) + 1), e3, marker='o')
plt.semilogy(range(1, len(e4) + 1), e4, marker='o')
plt.show()'''

# 3 mètode de la secant
'''res1, e1 = secant(-1, niter) # sol = -0,69
print(res1, f2(res1))

res2, e2 = secant(2, niter) # sol = 1,23
print(res2, f2(res2))

res3, e3 = secant(4, niter) # sol = 3,46
print(res3, f2(res3))

res4, e4 = secant(4, niter) # sol = -0,69
print(res4, f2(res4))

plt.semilogy(range(1, len(e1) + 1), e1, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e2) + 1), e2, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e3) + 1), e3, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e4) + 1), e4, marker='o') # convergència quadràtica
plt.show()
'''
# 4 aproximació de la derivada, com més petita la h convergeix més de pressa,
'''h = 1e-5
res1, e1 = newton_derivada(-1, niter, h) # sol = -0,69
print(res1, f2(res1))

res2, e2 = newton_derivada(2, niter, h) # sol = 1,23
print(res2, f2(res2))

res3, e3 = newton_derivada(3, niter, h) # sol = 3,46
print(res3, f2(res3))

res4, e4 = newton_derivada(2.5, niter, h) # sol = -0,69
print(res4, f2(res4))

plt.semilogy(range(1, len(e1) + 1), e1, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e2) + 1), e2, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e3) + 1), e3, marker='o') # convergència quadràtica
plt.semilogy(range(1, len(e4) + 1), e4, marker='o') # convergència quadràtica
plt.show()'''
