import numpy as np
import math as m
import matplotlib.pyplot as plt

def delay(arr, num, fill_value):
    result = np.empty_like(arr)
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result

#x is input
def diff_filter(x,T,L):
    num_coef = 2*L+1
    h_d = np.zeros(2*L+1)
    j = 0
    for n in range(-L,L+1):
        if n == 0:
            h_d[j] = 0
        else:
            h_d[j] = (1/T)*(((-1)**n)/n)
        j += 1
    blk = np.blackman(num_coef)
    h_blk = blk*h_d
    out = np.convolve(x,h_blk)
    return out

if __name__ == "__main__":
    print("Start Process...")
    t = np.arange(0,10,1/100) #1000 points from 0 to 10
    T = 1
    Fo = 0.15/T
    input = np.sin(2*np.pi*Fo*t)
    cos = np.cos(2*np.pi*Fo*t)
    
#differentiate and delay by num_coef
    out = diff_filter(input,100,5)
    L = 5
    num_coef = 2*L+1
#plot 
    plt.plot(t,input,label='original')
    plt.plot(t,(num_coef*1000)*out[num_coef-1:],color='red',label='differential') #mult by 1000 since t is 1000 length
    plt.legend(loc='upper right');plt.xlabel('Samples');plt.ylabel('Amplitude');plt.title('Differentiated Signal')
    plt.grid()
    
