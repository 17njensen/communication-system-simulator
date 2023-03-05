'''
Nik Jensen
'''
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

def generate_Ks(BnT, damp):
    theta_n = (BnT/(damp + 1/(4*damp)))
    k1_top = 4*damp*theta_n
    k2_top = 4*theta_n**2
    bottom = 1 + 2*damp*theta_n + theta_n**2
    k1 = k1_top/bottom
    k2 = k2_top/bottom
    return k1, k2

def err_check(prompt, prompt_answer):
    if type(prompt) != type(prompt_answer):
        print("INCORRECT ANSWER PROVIDED")
        sys.exit("Error message")
    else:
        return
if __name__ == "__main__":
    k1,k2 = generate_Ks(0.05, 1)
    print("C.2.1 Results:")
    print("k1 = ", k1)
    print("k2 = ", k2)

    k1,k2 = generate_Ks(0.1, 1.39)
    #k1 = 2*damp*wn
    #k2 = wn^2
    print("C.2.3 Results:")
    print("k1 = ", k1)
    print("k2 = ", k2)

    print("\nChoose Figure")
    print("1: C.2.4\n2: C.2.6\n3: C.2.8\n")
    choice = input("Choice: ")
    err_check(int(choice), 10) #error checks input 
    if (int(choice) > 3) or (int(choice) < 1):
        sys.exit("ERROR: Not a choice")
    input_cos = []
    output_cos = []
    freq = 2*m.pi/10
    phase = m.pi
    BnT = 0.05
    damp = 1
    conj_est = 0
    length = 200
    for n in range(0,length):
        n = n/2
        # input = m.exp(freq*n + phase)
        theta_input = freq*n + phase
        i_real = np.cos(theta_input)
        i_imag = np.sin(theta_input)

        input_cos.append(i_real)
        #phase detector
        if(n == 0):
            # p = input * conj_est
            #for C.2.4 and C.2.8
            if(choice == '1') or (choice == '2'):
                p_real = np.cos(theta_input + conj_est)
                p_imag = np.sin(theta_input + conj_est)
            if(choice == '1'):
                #C.2.4
                title = "Example C.2.4"
                arg = m.atan2(p_imag, p_real)
            elif(choice == '2'):
                #C.2.6
                title = "Example C.2.6"
                arg = p_imag
            elif(choice == '3'):
                #C.2.8
                title = "Example C.2.6"
                p_real = np.cos(theta_input - conj_est)
                p_imag = np.sin(theta_input - conj_est)
                arg = theta_input - conj_est

            k1,k2 = generate_Ks(BnT, damp)
            k1_arg = arg * k1
            k2_arg = arg * k2
            k2_arg_delay = k2_arg
            v_n_delay = 0
        else:
            #phase detector
            #for C.2.4 and C.2.8
            if(choice == '1') or (choice == '2'):
                p_real = np.cos(theta_input + conj_est)
                p_imag = np.sin(theta_input + conj_est)
            if(choice == '1'):
                #C.2.4
                arg = m.atan2(p_imag, p_real)
            elif(choice == '2'):
                #C.2.6
                arg = p_imag
            elif(choice == '3'):
                #C.2.8
                p_real = np.cos(theta_input - conj_est)
                p_imag = np.sin(theta_input - conj_est)
                arg = theta_input - conj_est

            #loop filter
            k1,k2 = generate_Ks(BnT, damp)
            k1_arg = arg * k1
            k2_arg = arg * k2
            k2_arg_delay = k2_arg + k2_arg_delay

        v_n = k1_arg + k2_arg_delay

        #DDS
        v_n_freq = v_n + freq + v_n_delay
        v_n_delay = v_n_freq

        # output = m.cos(v_n_delay) + j * m.sin(v_n_delay)
        output_real = m.cos(v_n_delay)
        output_imag = m.sin(v_n_delay)
        #C.2.4 and C.2.6
        if(choice == '1') or (choice == '2'):
            conj_est = -1*v_n_delay
        else:
            #C.2.8
            conj_est = v_n_delay
        output_cos.append(output_real)
    # plt.plot(output, 'red')
    x = np.arange(0,200)
    plt.plot(x/2,input_cos,color='blue',linestyle='--')
    plt.plot(x/2,output_cos, color='red')
    plt.title(title); plt.xlabel("Sample number"); plt.ylabel("Real part of sinusoids")
    plt.show()
