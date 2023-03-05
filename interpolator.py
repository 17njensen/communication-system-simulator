'''
Various Interpolators used to interpolate sinusoidal waves
'''
import numpy as np
import matplotlib.pyplot as plt
def linear_interpolate(x,T,x_prev):
    x_out = (x*T) + x_prev
    return x_out

#provide x[n+2] and mu, returns interpolated x[n]
def cubic_interpolate(x,mu):
    n = 0
    v_3 = mu*((x[n]*(1/6)) + (x[n-1]*(-1/2)) + (x[n-2]*(1/2)) + (x[n-3]*(-1/6)))
    v_2 = mu*(v_3 + ((x[n-1]*(-1/2)) + (x[n-2]*(-1)) + (x[n-3]*(1/2))))
    v_1 = mu*(v_2 + ((x[n]*(-1/6)) + (x[n-1]) + (x[n-2]*(-1/2)) + (x[n-3]*(-1/3))))
    x_out = x[n-2] + v_1
    return x_out
    
if __name__ == "__main__":
    print("Start Interpolation...")
    T = 5
    T0 = 1
    Fo = 0.15/T0
#create input
    t = np.arange(0,50,1/10)
    x = np.sin(2*np.pi*Fo*t)
    x_n = []
    for i in range(0,len(x),T):
        x_n.append(x[i])
#Linear Interpolate
    x_lin = []
    for i in range(1, len(x_n)):
        x_lin.append(linear_interpolate(x_n[i], T, x_n[i-1])/T)
#Cubic Interpolate
    mu = 0
    x_out = []
    for i in range(0+3,len(x_n)+1):
        x_out.append(cubic_interpolate(x_n[i-3:i],mu))
    t1 = np.linspace(0,50-1,len(x_out))
    t2 = np.linspace(0,50-0.5,len(x_n))
    t3 = np.linspace(0,50-0.5,len(x_lin))
#plot
    plt.figure()
    plt.plot(t,x,color='red',linewidth=2,label='original')
    plt.plot(t3,x_lin,color='black',label='linear interpolated')
    plt.plot(t2,x_n,'.',color='blue')
    plt.legend(loc='upper right')
    plt.title('Linear Interpolation Result')
    plt.xlabel('Amplitude');plt.ylabel('Samples')
    plt.grid()
    plt.figure()
    plt.plot(t,x,color='red',linewidth=2,label='original')
    plt.plot(t2,x_n,'.',color='blue')
    plt.plot(t1,x_out,color='black',label='cubic interpolated')
    plt.legend(loc='upper right')
    plt.title('Cubic Interpolation Result')
    plt.xlabel('Amplitude');plt.ylabel('Samples')
    plt.grid()
    plt.show()