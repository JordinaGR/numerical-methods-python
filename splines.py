import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy

dataX = np.array([0, 1, 3, 4, 5, 7])
dataY = np.array([1, 1.25, 1, 0.5, 0, 0.4])

def calculaDerivada(x,y):
    n = len(dataX) -1
    der = np.zeros(n+1)
    der[0] = (dataY[1] - dataY[0]) / (dataX[1] - dataX[0])
    i = np.arange(1, n)
    der[i] = (dataY[i+1] - dataY[i-1]) / (dataX[i+1] - dataX[i-1])
    der[n] = (dataY[-1] - dataY[-2]) / (dataX[-1] - dataX[-2])
    return der

def sistema(dataX, dataY):
    n = len(dataX)
    h = np.zeros(n-1)
    t = np.zeros(n-1)

    for i in range(n-1):
        h[i] = dataX[i+1]-dataX[i]
        t[i] = dataY[i+1]-dataY[i]

    lamb = np.zeros(n-2)
    mu = np.zeros(n-2)
    e = np.zeros(n-2)
    for i in range(len(mu)):
        lamb[i] = h[i+1] / (h[i+1] + h[i])
        mu[i] = h[i]/(h[i+1]+h[i])
        e[i] = 3*(h[i+1]*t[i]/h[i] + h[i]*t[i+1]/h[i+1])/(h[i+1] + h[i])
    A = np.zeros((n, n))

    for i in range(n-2):
        A[i, i+1] = 2
        A[i, i] = lamb[i]
        A[i, i+2] = mu[i]

    A[n-2, 0] = 2
    A[n-2, 1] = 1
    e = np.append(e, 3*t[0]/h[0])
    A[n-1, n-2] = 1
    A[n-1, n-1] = 2
    e = np.append(e, 3*t[n-2]/h[n-2])
    
    return np.linalg.solve(A, e)


# 2.
#dataY = np.zeros(len(dataX))
#der = np.zeros(len(dataX))
#dataY[2] = 1
#der[3] = 1
# --------------------------------------

#der = calculaDerivada(dataX, dataY) # derivada enunciat
der = sistema(dataX, dataY) # derivada per fer un spline cubic c2 amb la funcion spline cubic c1
splineCubicC1 = scipy.interpolate.CubicHermiteSpline(dataX, dataY, der)
# print('Coeficients del spline: ', splineCubicC1.c)

# spline cubic C2 calculat amb scipy
der2 = calculaDerivada(dataX, dataY)
splineCubicC2 = scipy.interpolate.CubicSpline(dataX, dataY) # funcio per calcular spline cubic c2

x = np.arange(0,7.02,0.05)
y1 = splineCubicC1(x)
dy1 = splineCubicC1(x, nu = 1)   # primera derivada del spline
d2y1 = splineCubicC1(x, nu = 2)  # segona derivada del spline
plt.plot(dataX, dataY,'k*')
plt.plot(x,y1,label = 'S')
plt.plot(x,dy1,label = 'dS')
plt.plot(x,d2y1, label = 'd2S')
plt.title('Spline cubic C1 i les seves derivades')
plt.grid()
plt.legend()
plt.show()

# copia pel segon gràfic
y2 = splineCubicC2(x)
dy2 = splineCubicC2(x, nu = 1)   # primera derivada del spline
d2y2 = splineCubicC2(x, nu = 2)  # segona derivada del spline
plt.plot(dataX, dataY,'k*')
plt.plot(x,y2,label = 'S')
plt.plot(x,dy2,label = 'dS')
plt.plot(x,d2y2, label = 'd2S')
plt.title('Spline cubic C2 i les seves derivades')
plt.grid()
plt.legend()
plt.show()



# 1. les derivades s'han aproximat fent derivació numèrica com ti/hi
# l'spline es continu i derivable C^1. la primera derivada és 
# continua però no derivable per tant la segona derivada és discontinua
# si, és l'esperada pq l'spline es C^1

# 2. la funcio de la base son tot zeros excepte un 1 a la posició i

# 3. hem resolt el sistema per trobar les derivades 

# 4. 
