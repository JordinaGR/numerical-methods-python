# -*- coding: utf-8 -*-
"""
Solution of the 1D parabolic equation u_t=nu*u_xx for x in (a,b), t in (0,T)
"""
import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt

def initialCondition(x):
    #return np.sin(np.pi*x)
    return 1-2*np.abs(x-0.5)+0.5*x


def explicitDirichlet(a,b,nu,T,ua,ub,u0Func,nOfIntervals,nOfSteps):
    At=T/nOfSteps
    Ax=(b-a)/nOfIntervals
    x=np.arange(a,b+Ax,Ax)
    u0=u0Func(x)
    nOfStepsPlot=np.round(nOfSteps/10)   
    U=np.zeros((nOfSteps+1,nOfIntervals+1))
    U[0,:]=u0
    U[:,0]=ua
    U[:,nOfIntervals]=ub
    r=nu*At/Ax**2
    print("Explicit method: r=%5.2f" % (r))
    plt.plot(x,U[0,:],'o-')
    for n in np.arange(0,nOfSteps):
        U[n+1,1:-1]=U[n,1:-1]+r*(U[n,0:-2]-2*U[n,1:-1]+U[n,2:])
        if np.mod(n+1,nOfStepsPlot)==0:
            plt.plot(x,U[n+1,:],'o-')
    plt.title("EXplicit method")
    plt.show()
    return x,U

def implicitDirichlet(a,b,nu,T,ua,ub,u0Func,nOfIntervals,nOfSteps):
    At=T/nOfSteps
    Ax=(b-a)/nOfIntervals
    x=np.arange(a,b+Ax,Ax)
    u0=u0Func(x)
    nOfStepsPlot=np.round(nOfSteps/10)   
    U=np.zeros((nOfSteps+1,nOfIntervals+1))
    U[0,:]=u0
    U[:,0]=ua
    U[:,nOfIntervals]=ub
    r=nu*At/Ax**2
    A=diags([-r, 1+2*r, -r], [-1, 0, 1], shape=(nOfIntervals-1,nOfIntervals-1)).toarray()
    print("Implicit method: r=%5.2f" % (r))
    plt.plot(x,U[0,:],'o-')
    for n in np.arange(0,nOfSteps):
        b=np.copy(U[n,1:-1])
        b[0]=b[0]+r*ua
        b[-1]=b[-1]+r*ub
        U[n+1,1:-1]=np.linalg.solve(A,b)
        if np.mod(n+1,nOfStepsPlot)==0:
            plt.plot(x,U[n+1,:],'o-')
    plt.title("IMplicit method")
    plt.show()
    return x,U

def explicitHomogeneousNeumann(a,b,nu,T,u0Func,nOfIntervals,nOfSteps):
    At=T/nOfSteps
    Ax=(b-a)/nOfIntervals
    x=np.arange(a,b+Ax,Ax)
    u0=u0Func(x)
    nOfStepsPlot=np.round(nOfSteps/10)   
    U=np.zeros((nOfSteps+1,nOfIntervals+1))
    U[0,:]=u0
    r=nu*At/Ax**2
    print("Explicit method: r=%5.2f" % (r))
    plt.plot(x,U[0,:],'o-')
    for n in np.arange(0,nOfSteps):
        U[n+1,0]=U[n,0]+r*(-2*U[n,0]+2*U[n,1])
        U[n+1,1:-1]=U[n,1:-1]+r*(U[n,0:-2]-2*U[n,1:-1]+U[n,2:])
        U[n+1,-1]=U[n,-1]+r*(2*U[n,-2]-2*U[n,-1])
        if np.mod(n+1,nOfStepsPlot)==0:
            plt.plot(x,U[n+1,:],'o-')
    plt.title("EXplicit method with homogeneous Neumann conditions")
    plt.show()
    return x,U

def implicitHomogeneousNeumann(a,b,nu,T,u0Func,nOfIntervals,nOfSteps):
    At=T/nOfSteps
    Ax=(b-a)/nOfIntervals
    x=np.arange(a,b+Ax,Ax)
    u0=u0Func(x)
    nOfStepsPlot=np.round(nOfSteps/10)   
    U=np.zeros((nOfSteps+1,nOfIntervals+1))
    U[0,:]=u0
    r=nu*At/Ax**2
    A=diags([-r, 1+2*r, -r], [-1, 0, 1], shape=(nOfIntervals+1,nOfIntervals+1)).toarray()
    A[0,1]=-2*r
    A[-1,-2]=-2*r
    print("Implicit method: r=%5.2f" % (r))
    plt.plot(x,U[0,:],'o-')
    for n in np.arange(0,nOfSteps):
        U[n+1,:]=np.linalg.solve(A,U[n,:])
        if np.mod(n+1,nOfStepsPlot)==0:
            plt.plot(x,U[n+1,:],'o-')
    plt.title("IMplicit method with homogeneous Neumann conditions")
    plt.show()
    return x,U


a=0
b=1
ua=initialCondition(a)
ub=initialCondition(b)
T=0.1
nu=1
nOfIntervals=10
nOfSteps=80
x,U = explicitDirichlet(a,b,nu,T,ua,ub,initialCondition,nOfIntervals,nOfSteps)
x,U = implicitDirichlet(a,b,nu,T,ua,ub,initialCondition,nOfIntervals,nOfSteps)
x,U = explicitHomogeneousNeumann(a,b,nu,T,initialCondition,nOfIntervals,nOfSteps)
x,U = implicitHomogeneousNeumann(a,b,nu,T,initialCondition,nOfIntervals,nOfSteps)











