import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial
from numpy.polynomial.polynomial import polymul, polyval
import scipy 
from scipy import integrate


# Probem settings
a = -1
b = 1
def f(x):
    return 1/(1+25*(x**2))


#------------------------------------------------------------
# Exemple: aproximem un conjunt de numPts punts (dataX, dataY)
# mitjançant un polinomi i amb el criteri de mínims quadrats
#------------------------------------------------------------
numPts = 101
degree = 2
dataX = np.linspace(-1,1, numPts)
dataY = f(dataX)

# Obtenim els coeficients del polinomi i l'avaluem en els punts que fem servir per dibuixar
'''poly = Polynomial.fit(dataX,dataY,degree) # minims quadrats
xx = np.linspace(-1,1, 101)
yy = poly(xx)

# Representem gràficament els resultats
plt.plot(dataX, dataY, 'k', label = 'f(x) = 1/(1+25x^2)')
plt.plot(xx, yy, label='grau = {}'.format(degree))
plt.title('Polinomi mínims quadrats')
plt.legend()
plt.show()
'''

# 1 fem minims quadrats 
'''
# ------------------------------------ calcul funcions i errors
degreeVec = [2, 4, 6, 8, 10]
plt.plot(dataX, dataY, 'k', label = 'f(x) = 1/(1+25x^2)')

for deg in degreeVec:
    poly = Polynomial.fit(dataX, dataY, deg)
    xx = np.linspace(-1,1, 101)
    yy = poly(xx)
    
    aproxY = poly(dataX)
    E = np.sqrt(np.sum( (dataY - aproxY)**2))
    func = lambda x: (f(x) - poly(x))**2
    r = np.sqrt(integrate.quad(func, a, b))
    print('\nMínims quadrats amb ' + str(numPts) + ' punts i un polinomi de grau ' + str(deg))
    print('E = ' + str(E) + ', int(res(x)) = ' + str(r[0]))
    
    plt.plot(xx, yy, label="grau = {}".format(deg))

plt.legend()
plt.show()
# no podem esperar el mateix resultat per qualsevol f(x)  ni qualsevol interval, per exemple f(x) = x
# ------------------------------------
'''

# 2 hem fet minims quadrats amb el nombre de  punts major al grau dona interpolació
# -----------------interpolacio------------------- calcul funcions i errors

'''degreeVec = [2, 4, 6, 8, 10]
plt.plot(dataX, dataY, 'k', label = 'f(x) = 1/(1+25x^2)')

for deg in degreeVec:
    dataX = np.linspace(-1,1, deg+1)
    dataY = f(dataX)

    poly = Polynomial.fit(dataX, dataY, deg)
    xx = np.linspace(-1,1, 101)
    yy = poly(xx)
    
    aproxY = poly(dataX)
    E = np.sqrt(np.sum( (dataY - aproxY)**2))
    func = lambda x: (f(x) - poly(x))**2
    r = np.sqrt(integrate.quad(func, a, b))
    print('\nMínims quadrats amb ' + str(numPts) + ' punts i un polinomi de grau ' + str(deg))
    print('E = ' + str(E) + ', int(res(x)) = ' + str(r[0]))
    
    plt.plot(xx, yy, label="grau = {}".format(deg))

plt.legend()
plt.show()
'''
# 1 punt mes que grau
# ------------------------------------


# p2 = np.polynomial.legendre.leg2poly([0,0,1])
# p3 = np.polynomial.legendre.leg2poly([0,0,0,1])
# print('Coef. del polinomi de Legendre de grau 2: ',p2)
# print('Coef. del polinomi de Legendre de grau 3: ',p3)
# # Calculem els coeficients del polinomi p2*p2
# prod = polymul(p2, p2)
# # i fem servir la funció Polynomial per convertir-lo en una funció polinòmica
# func = lambda x: Polynomial(prod)(x)
# # La funció scipy.integrate.quad aproxima integrals definides
# # (el resultat és una tupla amb el valor de la integral i una estimació de l'error)
# res = integrate.quad(func, a, b)
# print('<p2,p2> = ', res[0])
# prod = polymul(p2, p3)
# func = lambda x: Polynomial(prod)(x)
# res = integrate.quad(func, a, b)
# print('<p2,p3> = ', res[0])


# p1 = np.polynomial.chebyshev.cheb2poly([0,1])
# p4 = np.polynomial.chebyshev.cheb2poly([0,0,0,1])
# print('Coef. del polinomi de Chebyshev de grau 1: ',p1)
# print('Coef. del polinomi de Chebyshev de grau 4: ',p4)
# prod = polymul(p1, p1)
# # En aquest cas multipliquem el producte de polinomis pel pes (1 / sqrt(1+x^2))
# # per poder definir correctament el producte escalar
# func = lambda x: Polynomial(prod)(x) * (1/np.sqrt((1 - x**2)))
# res = scipy.integrate.quad(func, a, b)
# print('<p1,p1> = ', res[0])
# prod = polymul(p1, p4)
# func = lambda x: Polynomial(prod)(x) * (1/np.sqrt((1 - x**2)))
# res = scipy.integrate.quad(func, a, b)
# print('<p1,p4> = ', res[0])
