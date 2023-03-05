import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.random import randint
from interpolator import cubic_interpolate, linear_interpolate
from signals import srrc
from differentiate import diff_filter

def find_nearest(array, value): 
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

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
    print("symbol_timing_recovery process starting...")
    #init variables
    N = 16
    #for srrc
    alpha = 1
    lp = 10 
    Ts = 100 
    #for diff filter
    T_d = 100 
    L = 5
    T0 = 1
    Fo = 0.15/T0
    
    K1 = 2.6*10**-2
    K2 = 6.9*10**-4
    W1 = 0
    MU1 = 0
    V21 = 0
    NCO1 = 0
    dxout = []
    dxout_q = []
    xout = []
    xout_q = []
    u_k = []
    e_out = []
    
    Nsym = 100
#sender
    s_i = np.zeros(Nsym)
    s_q = np.zeros(Nsym)
    for n in range(0,Nsym):
        bit = randint(4)#QPSK
        s_1, s_2, LUT = PSK(4,bit) #from homework 4
        print("LUT = ",LUT)
        s_i[n] = s_1
        s_q[n] = s_2
    i_up = np.zeros((N*len(s_i),1))
    i_up[range(0,N*Nsym,N)] = s_i.reshape(Nsym,1)
    q_up = np.zeros((N*len(s_q),1))
    q_up[range(0,N*Nsym,N)] = s_q.reshape(Nsym,1)
    plt.figure()
    plt.plot(s_i)
    sig = srrc(alpha,N,lp,Ts)
    x_l = np.convolve(i_up[:,0],sig)
    x_q = np.convolve(q_up[:,0],sig)
    x = x_l
    
#receiver    
#run matched filter, return x_n and dx_n   
    print("size of original sig = ",len(sig))
    x_mf = np.convolve(x,sig)
    psrrc = np.sum(sig * sig) #power srrc
    x_n_not_delay = x_mf/psrrc #normalize 
    x_mfq = np.convolve(x_q,sig)
    x_n_not_delay_q = x_mfq/psrrc #normalize 
    
    dx_mf = np.convolve(x,sig)
    psrrc = np.sum(sig * sig) #power srrc
    dpx_n = dx_mf/psrrc #normalize 
    dx_mfq = np.convolve(x_q,sig)
    psrrc = np.sum(sig * sig) #power srrc
    dpx_nq = dx_mfq/psrrc #normalize 
    
    delta = np.zeros(2*L+1) #delay function
    delta[L] = 1
    
#differentiate, matched filter outputs
    x_n = np.convolve(x_n_not_delay,delta)
    x_n_q = np.convolve(x_n_not_delay_q,delta)

    dx_n = diff_filter(dpx_n,T_d,L)
    dx_n_q = diff_filter(dpx_nq,T_d,L)
#symbol timing recovery
    j = 3
    l = 1
    for n in range(0,(N-1)*len(x_n)):
        temp = NCO1 - W1
        if temp < 0:
            strobe = 1
            mu = NCO1/w
            nco = 1 + temp
        else:
            strobe = 0
            mu = MU1
            nco = temp
        if strobe == 0:
            e = 0
        else:
        #interpolate
            xi = cubic_interpolate(x_n[j-3:j],mu)       #INPHASE
            dxi = cubic_interpolate(dx_n[j-3:j],mu)
            xi_q = cubic_interpolate(x_n_q[j-3:j],mu)   #QUAD
            dxi_q = cubic_interpolate(dx_n_q[j-3:j],mu)
            # xi = linear_interpolate(x_n[l],N,(x_n[l-1])/N) #INPHASE
            # dxi = linear_interpolate(dx_n[l],N,dx_n[l-1]/N)
            # xi_q = linear_interpolate(x_n_q[l],N,(x_n_q[l-1])/N) #QUAD
            # dxi_q = linear_interpolate(dx_n_q[l],N,dx_n_q[l-1]/N)
        #compute error signal
            e = np.sign(xi) * dxi
            e_out.append(e)
        #a_k
            xout.append(xi)     #INPHASE
            dxout.append(dxi)
            xout_q.append(xi_q) #QUAD
            dxout_q.append(dxi_q)
            j += 1
            l += 1
    #loop filter
        v1 = K1*e
        v2 = V21 + K2*e
        v = v1 + v2
        w = v + 1/N
        
        V21 = v2
        NCO1 = nco
        MU1 = mu
        u_k.append(mu)
        W1 = w

    a = []
    b = []
    a_q = []
    b_q = []
    x_i_out = []
    x_q_out = []
#decisions 
    for i in range(0,len(xout),N):
        x_i_out.append(230264317*xout[i])#cubic
        x_q_out.append(230264317*xout_q[i])
        # x_i_out.append(xout[i]/(N-2))#linear
        # x_q_out.append(xout_q[i]/(N-2))
    for i in range(0,len(xout)):
        a.append(find_nearest(LUT.real, xout[i]))
        b.append(find_nearest(LUT.real, x_n[i]))
        a_q.append(find_nearest(LUT.imag, xout_q[i]))
        b_q.append(find_nearest(LUT.imag, dxout_q[i]))
#PLOT!!!
    plt.figure()
    plt.plot(x_i_out[2:],x_q_out[2:],'.')
    plt.plot(a,a_q,'x')
    plt.title('Scatter Plot for Cubic Interpolation')
    plt.figure()
    plt.plot(u_k);plt.title('Fractional Interval of Cubic Interpolation')
    plt.xlabel('Symbol Index');plt.ylabel('Fractional Interval u(k)')
    plt.figure()
    plt.plot(e_out);plt.title('Timing Error Signal for Cubic Interpolation')
    plt.xlabel('t/Ts');plt.ylabel('Timing Error e(t)')
    a.insert(0,0)
    a.insert(0,0)
    plt.figure()
    plt.plot(b[:],color='red',label='original')
    plt.plot(a[:],color='blue',label='interpolated')
    plt.legend(loc='upper right')
    plt.title('Decision Bits for Cubic Interpolation');plt.xlabel('Samples');plt.ylabel('Constellation Amplitude')
    plt.figure()
    plt.plot(b[:],color='red',label='original')
    plt.plot(a[:],color='blue',label='interpolated')
    plt.legend(loc='upper right')
    plt.title('Decision Bits for Cubic Interpolation (Zoomed in)');plt.xlabel('Samples');plt.ylabel('Constellation Amplitude')
    plt.xlim(0,250)
    plt.show()