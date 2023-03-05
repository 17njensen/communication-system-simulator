'''
Nik Jensen
'''
import numpy as np
import math as m
import matplotlib.pyplot as plt

def srrc(beta, N, lp, Ts):
    z = np.arange(-lp*N,(lp*N)+1)

    h = np.zeros(len(z))

    for i in range(len(z)):
        n = z[i]

        if (n == 0):
            h[i] = (1/np.sqrt(N))*((1-beta) + (4*beta/np.pi))
        elif (np.abs(n*4*beta) == N):
            h[i] = (beta/np.sqrt(2*N))*( (1+(2/np.pi))*np.sin(np.pi/(4*beta)) + (1-(2/np.pi))*np.cos(np.pi/(4*beta)) )
        else:
            h[i] = (1/np.sqrt(N))*( (np.sin(np.pi*n*(1-beta)/N)) + (4*beta*n/N)*(np.cos(np.pi*n*(1+beta)/N)) ) / ( (np.pi*n/N) * (1 - (4*beta*n/N)**2) )

    h = h/np.sqrt(N)

    return h

def nrz(N, lp, Ts):
    z = np.arange(-lp*N, lp*N+1)
    h = np.zeros(len(z))
    for i in range(0,Ts): #from 0 to Ts
        h[i] = 1/m.sqrt(Ts)
    return h

def man(N,lp,Ts):
    z = np.arange(-lp*N,lp*N+1)
    h = np.zeros(len(z))
    for i in range(0,int(Ts/2)):
        h[i] = 1/m.sqrt(Ts)
    for i in range(int(Ts/2),Ts):
        h[i] = -1/m.sqrt(Ts)
    return h