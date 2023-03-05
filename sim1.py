# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 15:51:28 2022

@author: Nik
"""

import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
import struct   # i
from symbol_timing_recovery import find_nearest
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
    fname = 'sim1_2022'
    # fname = 'sim3_2022'
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
    # twoints1 = fid.read(struct.calcsize(packints))
    # twoints = struct.unpack(packints,twoints1)
    # rows = twoints[0];
    # cols = twoints[1];
    # rows = 230#200
    # cols = 300
    nsym = 69000
    size = nsym * N
    # print("rows = ",rows)
    # print("cols = ",cols)
    image = fid.read().splitlines()
    fid.close()
    # image = np.frombuffer(image,dtype='uint8')
    x = ([float(z) for z in image])
    # x = image.reshape(cols,rows).T
    # plt.figure(2)
    # plt.imshow(255-x,cmap=plt.get_cmap('Greys'))
    # plt.title('Original image')
    # plt.show()
    
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
    x = x_r[range(0,len(x_r),N)]
    y = y_r[range(0,len(y_r),N)]
    residual = round((len(srrc)-1)/N)

#decison
    a_0 = []
    a_1 = []
    a_x = []
    a_y = []
    for i in range(0,len(x)):
        a_x.append(find_nearest(LUT[:,0], x[i]))
        a_y.append(find_nearest(LUT[:,1], y[i]))
        a_0.append(round(find_nearest(LUT[:,0], x[i])*math.sqrt(Enorm)))
        a_1.append(round(find_nearest(LUT[:,1], y[i])*math.sqrt(Enorm)))
#create 256 QAM bitmap
    a = np.zeros(shape=(M,2))
    a = create_bitmap(M)
    aout = []
    aout_no_UW = []
    # UW = []
    cols = 0
    loc = 0
    flag = 0 #this will get flagged if UW is found in the second row
    for i in range(0,len(a_0)):
    # for i in range(24,nsym+24):
        arr = np.where((a == (a_0[i],a_1[i])).all(axis=1))
        aout.append(arr[0][0])
    uw = find_unique_word(aout,residual)
    uwlen = len(uw)
    for i in range(0,len(aout)):
        #if UW found, previous data is one row? 
        #number of times this occurs is num of columns
        aout_no_UW.append(aout[i])
        if((i + 2) < len(aout)):
            # if(aout[i] == 162) & (aout[i+1] == 29) & (aout[i+2] == 92): #UW found
            if(aout[i] == uw[0]) & (aout[i+1] == uw[1]) & (aout[i+2] == uw[2]):
                if flag == 0:
                    UW = aout[i:i+uwlen]
                if flag == 1:
                    res = loc
                    rows = i - loc
                    res = res - rows - uwlen #find residual data
                loc = i + uwlen
                if (loc < len(aout)):
                    i = loc
                flag = flag + 1
                cols = cols + 1
    cols = cols - 1
    nsym = rows*cols 
    nsym_UW = (rows*cols) + (uwlen*cols)
    if (len(aout_no_UW) > (nsym)):
        #residual = len(aout_no_UW) - (row*col)
        #remove residual/2 from front
        out_image = aout_no_UW[res:nsym+res]
#decision output
    aout = np.reshape(aout[res:nsym_UW+res],(cols,rows+uwlen)).T
    out_image = np.reshape(out_image,(cols,rows)).T
    plt.figure()
    plt.plot(LUT[:,0],LUT[:,1],'.')
    plt.plot(a_x,a_y,'x')
    plt.title('Constellation results')
    plt.figure()
    plt.imshow(255-aout,cmap=plt.get_cmap('Greys'))
    plt.title('sim1_2022 Results')
    # plt.imshow(255-out_image,cmap=plt.get_cmap('Greys'))
    plt.show()
    
    
    