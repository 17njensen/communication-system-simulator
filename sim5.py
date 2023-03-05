# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 22:03:01 2022

@author: Nik
"""

import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
import struct   # i
from symbol_timing_recovery import find_nearest
from phase_recovery import generate_Ks
from interpolator import cubic_interpolate, linear_interpolate
from differentiate import diff_filter
from signals import srrc
packints = 'ii'

# Set parameters
f_eps        = 0.0 # Carrier frequency offset percentage (0 = no offset)
Phase_Offset = 0.0 # Carrier phase offset (0 = no offset)
t_eps        = 0.0 # Clock freqency offset percentage (0 = no offset)
T_offset     = 0.0 # Clock phase offset (0 = no offset)
Ts = 1             # Symbol period
N = 4              # Samples per symbol period

def find_unique_word(arr,end_length):
    uw = []
    out = []
    for i in range(len(arr)-1-end_length,len(arr)-1-300,-1):
        uw.append(arr[i])
    uw.reverse()
    # print(uw)
    flag = 0
    for i in range(0,len(uw)):
        for j in range(0,500):
            if (flag == 0) & (i+4 < len(uw)) & (j+4 < 500):
                if (uw[i] == arr[j]) & (uw[i+1] == arr[j+1]) & (uw[i+2] == arr[j+2]) & (uw[i+3] == arr[j+3]) & (uw[i+4] == arr[j+4]):
                    print(arr[j],arr[j+1],arr[j+2],arr[j+3],arr[j+4])
                    out.append(arr[j])
                    out.append(arr[j+1])
                    out.append(arr[j+2])
                    out.append(arr[j+3])
                    out.append(arr[j+4])
                    if (arr[j] == out[1]) & (arr[j+1] == out[2]) & (len(out) > 5):
                        flag = 1
                    i = i+5
                    j = j+5
                    
    print(out[:(len(out)-5)])
    return out[:(len(out)-5)]         
    
def create_bitmap(M):
    a = np.zeros(shape=(M,2))
    '''
    a = [[-15,-15],    #0
         [-15,-14],    #1
         [-15,-13],    #2
         ...
         ]
    '''
    k = 0
    for i in range(-15,15+1,2):
        for j in range(-15,15+1,2):
            a[k] = [i,j]
            k = k + 1
    # print(a)
    return a

if __name__ == "__main__":
    # fname = 'test_2022'
    # fname = 'sim1_2022'
    # fname = 'sim2_2022'
    fname = 'data/sim5_2022'
# Select modulation type
    # Use 8 bits per symbol and 256 square QAM
    B = 8;            # Bits per symbol (B should be even: 8, 6, 4, 2)
    # B = 4;
    bits2index = 2**np.arange(B-1,-1,-1)
    M = 2 ** B       # Number of symbols in the constellation
    Mroot = math.floor(2**(B/2))
    a = np.reshape(np.arange(-Mroot+1,Mroot,2),(2*B,1))
    b = np.ones((Mroot,1))
    LUT = np.hstack((np.kron(a,b), np.kron(b,a)))
    Enorm = np.sum(LUT ** 2) / M;
    LUT = LUT/math.sqrt(Enorm);
    
#Load in image
    fid= open(fname,'rb')
    nsym = 69000
    size = nsym * N
    image = fid.read().splitlines()
    fid.close()
    x = ([float(z) for z in image])
    
#received modulated signal-----------------------------------------------------

#Form I and Q
    Omega0 = math.pi/2*(1+f_eps)
    n = np.arange(len(x))
    C =  math.sqrt(2)*np.cos(Omega0*n + Phase_Offset)
    S = -math.sqrt(2)*np.sin(Omega0*n + Phase_Offset)
    nstd = 0   # for this test, there is no noise
    I = x * C
    Q = x * S
    
#matched filter
    EB = 0.7;  # Excess bandwidth
    To = (1+t_eps)
    Lp = 12;
    t = np.arange(-Lp*N,Lp*N+1) /N + 1e-8;  # +1e-8 to avoid divide by zero
    tt = t + T_offset;
    srrc = ((np.sin(math.pi*(1-EB)*tt)+ 4*EB*tt * np.cos(math.pi*(1+EB)*tt))
        /((math.pi*tt)*(1-(4*EB*tt)**2)))
    srrc = srrc/math.sqrt(N);
#TED stuff
    psrrc = np.sum(srrc * srrc) #power srrc
    x_r = np.convolve(srrc,I) #inphase
    x_r = x_r/psrrc
    y_r = np.convolve(srrc,Q) #quad
    y_r = y_r/psrrc
    
    dx_r = np.convolve(srrc,I) #inphase
    dx_r = dx_r/psrrc
    dy_r = np.convolve(srrc,Q) #quad
    dy_r = dy_r/psrrc
#differentiate
    L = 5
    delta = np.zeros(2*L+1) #delay function - needs to change?
    delta[L] = 1
    x_r = np.convolve(x_r,delta) #inphase
    y_r = np.convolve(y_r,delta) #inphase
    
    dx_r = diff_filter(dx_r,100,5) #diff inphase
    dy_r = diff_filter(dy_r,100,5) #diff quad

    #perform the roation, using downsampled x and y
    # plt.figure()
    # plt.plot(x,y,'.')
    x_k = x_r
    y_k = y_r
    dx_k = dx_r
    dy_k = dy_r

    residual = round((len(srrc)-1)/N)
    #check the rotations for all 4 angles
    
#PED
    ep = 0
    theta = 0
    e_k2_prev = 0
    theta_out = []
    theta_out.append(theta)
    a_0 = []
    a_1 = []
#generate Ks
    BnT = 0.0225
    damp = 1/(np.sqrt(2))
    k1,k2 = generate_Ks(BnT,damp)
    # k1 = 2.6*10**-2
    # k2 = 6.9*10**-4
    k0 = 1
    omega = 0
    v_sum_prev = 0
#TED
    et = 0
    W1 = 0
    MU1 = 0
    V21 = 0
    NCO1 = 0
    #outputs
    eout = []
    xout = []
    dxout = []
    yout = []
    dyout = []
    u_k = []
    wout = []
    a_x = []
    a_y = []
    k = 1
    FX1 = 0; FX2 = 0; FX3 = 0; FX4 = 0; FX5 = 0;
    FDX1 = 0; FDX2 = 0; FDX3 = 0; FDX4 = 0; FDX5 = 0;
    FY1 = 0; FY2 = 0; FY3 = 0; FY4 = 0; FY5 = 0;
    FDY1 = 0; FDY2 = 0; FDY3 = 0; FDY4 = 0; FDY5 = 0;
    tempFx = 0; tempFdx = 0;
    tempFy = 0; tempFdy = 0
    jj = 3
    l = 1
    for n in range(0,len(x_k)):#(N)*(len(x_k)-3)):
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
            et = 0
        else:
        #interpolate
            tempFx = -0.5*x_k[n]
            tempFdx = -0.5*dx_k[n]
            tempFy = -0.5*y_k[n]
            tempFdy = -0.5*dy_k[n]
            
            VF2 = -tempFx + FX1 + FX2 - FX3;
            VF1 = tempFx + (- FX1 + FX4) + FX2 + FX3;
            VF0 = FX5;
            xi = (VF2*mu + VF1)*mu + VF0; # Interpolated signal
           
            VF2 = -tempFdx + FDX1 + FDX2 - FDX3;
            VF1 = tempFdx + (- FDX1 + FDX4) + FDX2 + FDX3;
            VF0 = FDX5;
            dxi = (VF2*mu + VF1)*mu + VF0; #% Interpolated signal
            
            VF2 = -tempFy + FY1 + FY2 - FY3;
            VF1 = tempFy + (- FY1 + FY4) + FY2 + FY3;
            VF0 = FY5;
            yi = (VF2*mu + VF1)*mu + VF0; # Interpolated signal
            
            VF2 = -tempFdy + FDY1 + FDY2 - FDY3;
            VF1 = tempFdy + (- FDY1 + FDY4) + FDY2 + FDY3;
            VF0 = FDY5;
            dyi = (VF2*mu + VF1)*mu + VF0;# % Interpolated signal
            
            
            a0 = find_nearest(LUT[:,0], xi) #decision for I
            a1 = find_nearest(LUT[:,1], yi) #decision for Q
            a_x.append(a0)
            a_y.append(a1)
            a_0.append(round(a0*math.sqrt(Enorm)))
            a_1.append(round(a1*math.sqrt(Enorm)))
            et = a0*dxi + a1*dyi
            eout.append(et)
            jj += 1
            l += 1
    #loop filter TED
        v1t = k1*et
        v2t = V21 + k2*et
        vt = v1t + v2t
        w = vt + 1/N
    #mod TED
        FX3 = FX2;
        FX2 = FX1;
        FX1 = tempFx;
        FX5 = FX4;
        FX4 = x_k[n];
        FDX3 = FDX2;
        FDX2 = FDX1;
        FDX1 = tempFdx;
        FDX5 = FDX4;
        FDX4 = dx_k[n];
        
        FY3 = FY2;
        FY2 = FY1;
        FY1 = tempFy;
        FY5 = FY4;
        FY4 = y_k[n];
        FDY3 = FDY2;
        FDY2 = FDY1;
        FDY1 = tempFdy;
        FDY5 = FDY4;
        FDY4 = dy_k[n];
        
        V21 = v2t
        NCO1 = nco
        MU1 = mu
        u_k.append(mu)
        W1 = w
#end

    a0_prev = [] 
    a1_prev = []    
    aout_prev = []
#find rows and columns
    a = np.zeros(shape=(M,2))
    a = create_bitmap(M)

    uwlen = 30
    flag = 0 #this will get flagged if UW is found in the second row
    loc = 0
    cols = 0
    j = 0
    #find dimensions - rows are always 230
    for i in range(residual,len(a_0)): #start (len of srrc)/4 after
        j = j + 1
        if(j % 230 == 0):
            j = 0
            cols = cols + 1
            
    rows = 200
    nsym = rows*cols 
    nsym_UW = (rows*cols) + (30*cols)
    # cols = 210
    aout = []
    for i in range(0,len(a_0)):
        arr = np.where((a == (a_0[i],a_1[i])).all(axis=1))
        aout.append(arr[0][0])
    aout = np.reshape(aout[residual-9:nsym_UW+residual-9],(cols,rows+uwlen)).T
    # aout = np.reshape(aout_prev[residual:nsym_UW+residual],(cols,rows+uwlen)).T
#decision output
    plt.figure()
    plt.plot(LUT[:,0],LUT[:,1],'.')
    plt.plot(a_x,a_y,'x')
    plt.title('Constellation results')
    plt.figure()
    plt.plot(u_k)
    plt.title('Sim5_2022 Fractional Interval'); plt.xlabel("Symbol Index");plt.ylabel("Fractional Interval u(k)")
    plt.figure()
    plt.plot(eout)
    plt.title("Sim5_2022 Timing Error Signal"); plt.xlabel("Symbols index");plt.ylabel("Error e(t)")
    plt.figure()
    plt.imshow(255-aout,cmap=plt.get_cmap('Greys'))
    plt.title("sim5_2022 Results")
    # plt.imshow(255-out_image,cmap=plt.get_cmap('Greys'))
    plt.show()
    
