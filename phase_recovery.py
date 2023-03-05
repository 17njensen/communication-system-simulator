'''
Nik Jensen
'''
from error_communication import PSK, find_nearest
from pll import generate_Ks
from signals import srrc
from numpy.random import randint, normal
import numpy as np
import matplotlib.pyplot as plt
import math
from Encoder_Decoder import encode_qpsk, encode_qpsk_table

def conv_SRRC(i_up, q_up, alpha, N, lp, Ts):
    sig = srrc(alpha,N,lp,Ts)
    print("size of original sig = ",len(sig))
    i_out = np.convolve(i_up,sig)
    q_out = np.convolve(q_up,sig)
    return i_out,q_out, sig

def generate_Ks(BnT, damp):
    theta_n = (BnT/(damp + 1/(4*damp)))
    k1_top = 4*damp*theta_n
    k2_top = 4*theta_n**2
    bottom = 1 + 2*damp*theta_n + theta_n**2
    k1 = k1_top/bottom
    k2 = k2_top/bottom
    return k1, k2

def points_to_bits_qpsk(a_0,a_1):
    if a_0 == 1 and a_1 == 1:
        return 0b00
    elif a_0 == -1 and a_1 == 1:
        return 0b01
    elif a_0 == -1 and a_1 == -1:
        return 0b10
    elif a_0 == 1 and a_1 == -1:
        return 0b11

def bits_to_points_qpsk(bits):
    if int(bits) == 0:
        return 1, 1
    elif int(bits) == 1:
        return -1, 1
    elif int(bits) == 2:
        return -1, -1
    elif int(bits) == 3:
        return 1, -1

def decode_qpsk(bits):
    # print("Decoding given data")
    out_bits = 0b00
    #if first bit == 1 (i.e. 1011 or 1000 or 1100), flip all bits
    if (bits >= 8):
        print("x = ",bin(bits))
        bits = (~bits & 0xF) #added 0xF due to ~bits flips to 2's complement
        print("~x = ",bin(bits))
    #if 0000 or 0101, 00
    if bits == 0b0000 or bits == 0b0101:
        out_bits = 0b01#00
    #if 0001 or 0111, 01
    if bits == 0b0001 or bits == 0b0111:
        out_bits = 0b10#01
    #if 0010 or 0100, 10
    if bits == 0b0010 or bits == 0b0100:
        out_bits = 0b11#10
    #if 0011 or 0110, 11
    if bits == 0b0011 or bits == 0b0110:
        out_bits = 0b00#11
    return out_bits

if __name__ == "__main__":
    print("Begin phase recovery system...")
    theta = 0
    N = 3
    omega = 0  # 10 #freq of modulators, omega = woT
    v_sum_prev = 0
    e_k2_prev = 0
    iterations = 300
    BnT = 0.02
    damp = 1/(np.sqrt(2))
    # k1,k2 = generate_Ks(BnT,damp)
    k1 = 2.6*10**-2
    k2 = 6.9*10**-4
    k0 = 1
    alpha = 1
    lp = 10
    Ts_srrc = 100
    Nsym = 300
    theta_init = 30*np.pi/180
    theta_out = []
    err_sig = []
    a_k = []
    dec_out = []
    theta_out.append(theta)
    err_sig.append(0)
#Generate random QPSK symbols 
    s_i = np.zeros(Nsym)
    s_q = np.zeros(Nsym)
    dk = []
    bit_plt1 = []
    bit_plt = []
    for n in range(0,Nsym):
        bit_in = randint(4)#QPSK
        # bit_plt1.append(bit_in)
        # print("original bit = ",bit_in)
        # #encode
        # if n == 0:
        #     #create random encoded bits 
        #     encode_bit = randint(4)#for encoder
        #     if encode_bit == bit_in: #for some reason it's not always random? 
        #         encode_bit = randint(4)
        #     print("random bit =",encode_bit)
        #     en_bit_b = list(bin(encode_bit).replace("0b",""))
        #     if encode_bit == 0 or encode_bit == 1:
        #         en_bit_b.insert(0,'0')
        # else:
        #     en_bit_b = []
        #     en_bit_b.append(dk[n-1])
        #     en_bit_b.append(dk[n]) 
        # #convert actual bits 
        # bit_b = list(bin(bit_in).replace("0b","")) #convert to binary, convert 0000 -> '0','0','0','0'
        # if bit_in == 0 or bit_in == 1:
        #     bit_b.insert(0,'0')#add zero to beginning if 1 or 0 to make 0x0 -> 0x00
        # print(bit_b[0] + bit_b[1] + en_bit_b[0] + en_bit_b[1])
        # encode_out = encode_qpsk(bit_b[0],bit_b[1],en_bit_b[0],en_bit_b[1])
        # en_list = list(encode_out)
        # dk.append(en_list[0])#d_2k-2 - n
        # dk.append(en_list[1])#d_2k-1 - n+1
        # bit = int(''.join(en_list[2:4]),2) #d_2k and d_2k+1 combined
        # bit_plt.append(bit)
        # print("encoded bit = " + str(bit) + "\n")
        # #encoded result returns b,b-1,d,d-1
        bit = bit_in
        s_1, s_2, LUT = PSK(4,bit) #from homework 4
        print("LUT = ",LUT)
        s_i[n] = s_1
        s_q[n] = s_2
    plt.plot(bit_plt)
    plt.plot(bit_plt1)
    plt.show()
#Upsample by N
    i_up = np.zeros((N*len(s_i),1))
    q_up = np.zeros((N*len(s_q),1))
    print(i_up.size)
    i_up[range(0,N*Nsym,N)] = s_i.reshape(Nsym,1)
    q_up[range(0,N*Nsym,N)] = s_q.reshape(Nsym,1)
#Pass through the SRRC transmitter pulse shaping filter on I and Q branches
    i_p, q_p, srrc_sig = conv_SRRC(i_up[:,0], q_up[:,0],alpha,N,lp,Ts_srrc)
    psrrc = np.sum(srrc_sig * srrc_sig)
    print('srrc power=',psrrc)
    # rotate the signal
    ip1 = np.cos(theta_init)*i_p - np.sin(theta_init)*q_p
    qp1 = np.sin(theta_init)*i_p + np.cos(theta_init)*q_p
    i_p = ip1
    q_p = qp1
#Modulate I by cos and Q by sin
    r_i = i_p;
    r_q = q_p;#redundant but I'm keeping it for now
#transmitted data has been completed
    print("len(r_i) = ",len(r_i))


#------------------------------------------------------------------------------------------------
#receiver processing

#matched filter on I/Q branches - due to symmetry of SRRC, this should(?) be the same as pulse shaping step since we shift to be causal
    x, y, srrc_sig = conv_SRRC(r_i, r_q,alpha,N,lp,Ts_srrc)
    x = x/psrrc #normalize 
    y = y/psrrc #normalize
#downsample
    x_k = x[round(len(srrc_sig)-1):len(x):N] #start at shift of srrc_sig
    y_k = y[round(len(srrc_sig)-1):len(x):N]
    # print("x_k = ", x_k)
    plt.plot(x_k,y_k,'.',color='blue')
    plt.plot(s_i,s_q,'x',color='red')
    plt.show()
    e = []
    a0_out = []
    a1_out = []
    print("length of x_k =",len(x_k))
    for i in range(0,len(x_k)):
    #compute the rotation parameters
        C = np.cos(theta)
        S = np.sin(theta)
    #perform the roation, using downsampled x and y
        xr = C*x_k[i] - S*y_k[i]
        yr = S*x_k[i] + C*y_k[i]
    #evaluate the ML PED and DECODE
        #these are our encoded bits from receiver
        a_0 = find_nearest(LUT.real, xr) #decision for I
        a_1 = find_nearest(LUT.imag, yr) #decision for Q
        a0_out.append(a_0)
        a1_out.append(a_1)
        # p2b = bin(points_to_bits_qpsk(round(a_0),round(a_1))).replace("0b","")
        # if int(p2b) == 0 or int(p2b) == 1:
        #     p2b = '0' + p2b
        # a_k.append(p2b) #link back to integer and return binary
        # a_k_list = list(a_k[i])
        # d_2k_2 = dk[i]#0 
        # d2_2k_1 = dk[i+1]#1 
        # d3_2k = a_k_list[0] #from a_k, #first bit of a_k[i]
        # d4_2k_1 = a_k_list[1] #from a_k #second bit of a_k[i]
        # str_dec = str(d_2k_2) + str(d2_2k_1) + str(d3_2k) + str(d4_2k_1)
        # est_dec = int(str_dec,2)
        # dec = decode_qpsk(est_dec)
        # dec_out.append(dec) #to be plotted

        # a_0,a_1 = bits_to_points_qpsk(dec)

    #find e_k
        e_k = yr*a_0 - xr*a_1
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
        
#plot
    plt.figure()
    plt.plot(a0_out,a1_out,'.')
    plt.title('aout')
    plt.show()
    plt.figure()
    plt.plot(dec_out[0:Nsym],color = 'blue',label="Decoder")
    plt.plot(bit_plt1,color='green',label="Original")
    plt.legend(loc="upper right")
    plt.xlabel("Symbols");plt.ylabel("Amplitude");plt.title("Decoder Results (Original Bits vs Decoded)")
    plt.figure()
    plt.plot(e,color='blue')
    plt.title("phase error signal")
    plt.figure()
    plt.plot(theta_out, color='red')
    plt.title("phase estimate")
    plt.show()
