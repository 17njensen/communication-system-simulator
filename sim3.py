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
    fname = 'data/sim3_2022'
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
    
#Load and display the image
    fid= open(fname,'rb')
    nsym = 69000
    size = nsym * N
    image = fid.read().splitlines()
    fid.close()
    x = ([float(z) for z in image])
    
#received modulated signal
    # x = x.flatten('F')
#Form I and Q
    Omega0 = math.pi/2*(1+f_eps)
    n = np.arange(len(x))
    C =  math.sqrt(2)*np.cos(Omega0*n + Phase_Offset)
    S = -math.sqrt(2)*np.sin(Omega0*n + Phase_Offset)
    nstd = 0   # for this test, there is no noise
    I = np.zeros(len(x))
    Q = np.zeros(len(x))
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
    x_r = np.convolve(srrc,I)
    y_r = np.convolve(srrc,Q)

#downsample
    x_k = x_r[range(0,len(x_r),N)]
    y_k = y_r[range(0,len(y_r),N)]
    residual = round((len(srrc)-1)/N)

#phase recovery
    e = []
    theta = 0
    e_k2_prev = 0
    theta_out = []
    theta_out.append(theta)
    a_0 = []
    a_1 = []
    #generate Ks
    BnT = 0.6
    damp = 1/(np.sqrt(2))
    k1,k2 = generate_Ks(BnT,damp)
    # k1 = 2.6*10**-2
    # k2 = 6.9*10**-4
    k0 = 1
    omega = 0
    v_sum_prev = 0
    a_x = []
    a_y = []
    for i in range(0,len(x_k)):
    #compute the rotation parameters
        C = np.cos(theta)
        S = np.sin(theta)
    #perform the roation, using downsampled x and y
        xr = C*x_k[i] - S*y_k[i]
        yr = S*x_k[i] + C*y_k[i]
    #evaluate the ML PED and DECODE
        #these are our encoded bits from receiver
        a0 = find_nearest(LUT[:,0], xr) #decision for I
        a1 = find_nearest(LUT[:,1], yr) #decision for Q
        a_x.append(find_nearest(LUT[:,0], xr))
        a_y.append(find_nearest(LUT[:,1], yr))
        a_0.append(round(a0*math.sqrt(Enorm)))
        a_1.append(round(a1*math.sqrt(Enorm)))
    #find e_k
        e_k = yr*a0 - xr*a1
        e.append(e_k)
    #apply the loop filter #refer to page 509 of pdf
        e_k1_ = e_k*k1
        e_k2 = e_k*k2
        e_k2_ = e_k2 + e_k2_prev
        e_k2_prev = e_k2_
        v = e_k1_ + e_k2_

    #DDS (baseband) #refer to Figure C.2.3
        v_k0 = v * k0 #assuming k0 = 1
        v_sum = v_k0 + omega + v_sum_prev #theta? previous v_sum
        v_sum_prev = v_sum
        theta = -v_sum_prev

        theta_out.append(theta*180/np.pi)
    print("size of theta_out =", len(theta_out))
    #end

#decison
    a0_prev = []
    a1_prev = []
    aout_prev = []
#find rows and columns (has some extra computations that arent necessary)
    a = np.zeros(shape=(M,2))
    a = create_bitmap(M)
    for i in range(0,len(x_k)):
        a0_prev.append(round(find_nearest(LUT[:,0], x_k[i])*math.sqrt(Enorm)))
        a1_prev.append(round(find_nearest(LUT[:,1], y_k[i])*math.sqrt(Enorm)))
    # for i in range(0,len(a_0)):
        arr = np.where((a == (a0_prev[i],a1_prev[i])).all(axis=1))
        aout_prev.append(arr[0][0])
    # uw = find_unique_word(aout_prev,residual)
    flag = 0 #this will get flagged if UW is found in the second row
    loc = 0
    cols = 0
    for i in range(residual,len(aout_prev)):
        if (i % 230 == 0): 
            cols = cols + 1
            uw = aout_prev[i-30:i]
    uwlen = len(uw)
    rows = 200
    cols = cols - 1
    nsym = rows*cols 
    nsym_UW = (rows*cols) + (uwlen*cols)
    
    aout = []
    for i in range(0,len(a_0)):
        arr = np.where((a == (a_0[i],a_1[i])).all(axis=1))
        aout.append(arr[0][0])
    
#create 256 QAM bitmap
    aout_no_UW = []
    # UW = []
    # for i in range(0,len(a_0)):
    # # for i in range(24,nsym+24):
    #     arr = np.where((a == (a_0[i],a_1[i])).all(axis=1))
    #     aout.append(arr[0][0])
    # uw = find_unique_word(aout,residual)
    # uwlen = len(uw)
    # if (len(aout_no_UW) > (nsym)):
    #     #residual = len(aout_no_UW) - (row*col)
    #     #remove residual/2 from front
    #     out_image = aout_no_UW[res:nsym+res]
#decision output
    aout = np.reshape(aout[residual:nsym_UW+residual],(cols,rows+uwlen)).T
    # aout = np.reshape(aout[residual:len(aout)-1],(cols,rows+uwlen)).T
    # out_image = np.reshape(out_image,(cols,rows)).T
    plt.figure()
    plt.plot(LUT[:,0],LUT[:,1],'.')
    plt.plot(a_x,a_y,'x')
    plt.title('Constellation results')
    plt.figure()
    plt.plot(e)
    plt.title("Sim3_2022 Error Signal"); plt.xlabel("Symbols index");plt.ylabel("Error")
    plt.figure()
    plt.plot(theta_out)
    plt.title("Sim3_2022 Phase Estimate"); plt.xlabel("Symbols index");plt.ylabel("Phase")
    plt.figure()
    plt.imshow(255-aout,cmap=plt.get_cmap('Greys'))
    plt.title("sim3_2022 Results")
    # plt.imshow(255-out_image,cmap=plt.get_cmap('Greys'))
    plt.show()
    
