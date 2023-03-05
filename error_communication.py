'''
Nik Jensen
'''
import numpy as np
from numpy.random import randint, normal
import matplotlib.pyplot as plt
import math
from scipy import special as sp

#Q-function
def qfunc(x):
    return 0.5*sp.erfc(x/np.sqrt(2))

#given value, find the nearest element in the given array. Returns value found
def find_nearest(array, value): 
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#decision for CCITT
def CCITT(bits_in):
    M = 4
    LUT = np.zeros(2*M,dtype = 'complex_')
    t = np.linspace(0,3*np.pi/2,M) #even spaced from 0:2pi over interval M
    t_1 = np.linspace(np.pi/4,7*np.pi/4,M)
    x_inner = np.cos(t_1)*(2/(np.sqrt(2)))
    y_inner = np.sin(t_1)*(2/(np.sqrt(2)))
    x_outer = 3*np.cos(t)
    y_outer = 3*np.sin(t)
    x = np.append(x_inner,x_outer)
    y = np.append(y_inner,y_outer)
    for i in range(0,len(x)):
        LUT[i] = x[i] + 1j*y[i]
    for i in range(0,len(x)):
        if bits_in == i:
            ampa_I = x[i]
            ampa_Q = y[i]
    return ampa_I,ampa_Q, LUT

#decision for M-PSK
def PSK(M,bits_in):
    LUT = np.zeros(M,dtype = 'complex_')
    b = (2*np.pi)/M
    t = np.linspace(0,b*(M-1),M) #even spaced from 0:2pi over interval M
    if M == 4:
        a = (2/(np.sqrt(2)))
    else:
        a = 1
    x = np.cos(t+np.pi/4)*a
    y = np.sin(t+np.pi/4)*a
    for i in range(0,M):
        LUT[i] = x[i] + 1j*y[i]
    # #assuming 4 bits of data
    for i in range(0,len(x)):
        if bits_in == i:
            ampa_I = x[i]
            ampa_Q = y[i]
    return ampa_I,ampa_Q, LUT

if __name__ == "__main__":
    Eb = 1
    No = 1
    LUT = np.array([-1,1]) #bpsk - M = 2
    M = 4
    Eavg = np.log2(M)*Eb
    maxsymbolerrorcount = 10
    Psymerrorarray = []
    Pbiterrorarray = []
    P = []
    startSNRdB = 0
    endSNRdB = 11
    SNRdBstep = 1
    isnr = 0
    bits_per_symb = int(input("bits per symbol = "))
    print("Start loop")
    print("symbolerrorcount = ", 0)
    print("maxsymbolerrorcount = ", maxsymbolerrorcount)
    for SNRdB in range(startSNRdB,endSNRdB,SNRdBstep):
        print(SNRdB)
        No = Eb/(10**(SNRdB/10)) 
        P.append((6/4)*qfunc(np.sqrt((2/3)*Eavg/No)) + (6/4)*qfunc(np.sqrt(2*Eavg/No)))
        # P.append(qfunc(2*np.sqrt(Eb/No)))
        sigma_sqr = No/2
        sigma = math.sqrt(sigma_sqr)
        symbolerrorcount = 0
        numsymbols = 0
        numbits = 0
        biterrorcount = 0
        while(symbolerrorcount < maxsymbolerrorcount):
            # print("looping")
            numsymbols += 1
            numbits += bits_per_symb #number of bits per symbol ([-1] = 0, [1] = 1 --> 1 bit)

            #transmitter
            bit = randint(2)#BPSK
            # bit = randint(4)#QPSK
            # bit = randint(8)#8PSK
            s = LUT[bit]  #for bpsk
            # s, s_q, LUT = PSK(4,bit) 
            # s, s_q, LUT = CCITT(bit)

            #channel w/ noise
            n = normal(0,sigma) #assuming number of bits, genreate "in each coordinate direction"
            n_q = normal(0,sigma)
            r = s + n
            # r_q = s_q + n_q

            #receiver
            # shat = find_nearest(LUT, r+1j*r_q)
            shat = find_nearest(LUT, r)#bpsk
            # shat_q = find_nearest(LUT, r_q)

            # if(shat.real != s) or (shat.imag != s_q):
            if(shat != s):#bpsk
                symbolerrorcount += 1
                print("\nSymbolerrorcount = ", symbolerrorcount)
                print("s = ",s)
                print("r = ",r)
                print("n = ",n)
                print("shat = ",shat)
                biterrorcount += 1 #the number of bits in error
        Psymerror = symbolerrorcount/numsymbols
        Pbiterror = biterrorcount/numbits
        Psymerrorarray.append(Psymerror)
        Pbiterrorarray.append(Pbiterror)
    #plots should be done on a logarithmic y axis
    plt.figure(); plt.clf()
    plt.semilogy(range(startSNRdB, endSNRdB, SNRdBstep), Pbiterrorarray, label="P(bit error)")
    plt.semilogy(range(startSNRdB, endSNRdB, SNRdBstep), Psymerrorarray, label="P(symbol error)")
    plt.semilogy(range(startSNRdB, endSNRdB, SNRdBstep), P, label="Predicted")
    plt.legend()
    plt.title("BPSK Error Simulation Results")
    plt.xlabel("Eb/No (dB)"); plt.ylabel("Probability of Error")
    plt.show()