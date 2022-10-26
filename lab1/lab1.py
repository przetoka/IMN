from cmath import sqrt
import matplotlib.pyplot as plt
import math
import numpy as np

analiticx = np.linspace(0.0, 5.0)
analiticy = np.exp(-1*analiticx)

def draw(x, y, filename,  desc, y_err = [],  err = 1, filename_err = "err.png",  additional = 0):
    plt.figure()
    plt.plot(analiticx, analiticy, color='yellow', linewidth = 2)
    plt.plot(x[0], list(y[0]), color='green', linestyle='dashed', linewidth = 2)
    plt.plot(x[1], list(y[1]),  color='red', linewidth = 1)
    plt.scatter(x[2], list(y[2]), color='blue', marker = '.', linewidth = 1)
    plt.legend(['analitic', 'dt = 0.01', 'dt = 0.1', 'dt = 1.0'])
    plt.ylabel('y(t)')
    plt.xlabel('t')
    plt.title(desc+ 'rozwiazanie')
    plt.savefig(filename)
    
    if err == 1:
        plt.figure()
        if additional == 1:
            plt.subplot(211)
        plt.plot(x[0], list(y_err[0]), color='green', linewidth = 1)
        plt.plot(x[1], list(y_err[1]),  color='red', linewidth = 1)
        plt.plot(x[2], list(y_err[2]), color='blue', linewidth = 1)
        plt.legend(['dt = 0.01', 'dt = 0.1', 'dt = 1.0'])
        plt.ylabel('y(t)')
        plt.xlabel('t')
        plt.title(desc+ 'error')
        if additional == 1:
            plt.subplot(212)
            plt.plot(x[0], list(y_err[0]), color='green', linewidth = 1)
            plt.plot(x[1], list(y_err[1]),  color='red', linewidth = 1)
            plt.legend(['dt = 0.01', 'dt = 0.1'])
            plt.ylabel('y(t)')
            plt.xlabel('t')
        plt.savefig(filename_err)

def euler_method(deltat, y0 = 1.0, lam = -1, tmax = 5.0, tmin  = 0):
    t = tmin
    y = y0
    result = []
    while t < tmax + deltat:
        result.append([y, y- np.exp(-1*t)])
        y += deltat * lam * y
        t += deltat
    return result

def rk2_method(deltat, y0 = 1.0, lam = -1, tmax = 5.0, tmin  = 0):
    t = tmin
    y = y0
    k1, k2 = 0, 0
    result = []
    while t <  tmax + deltat:
        result.append([y, y- np.exp(-1*t)])
        k1 = lam* y
        k2 = lam* (y + deltat* k1)
        y += deltat * (k1 + k2)/2
        t += deltat
    return result

def rk4_method(deltat, y0 = 1.0, lam = -1, tmax = 5.0, tmin  = 0):
    t = tmin
    y = y0
    k1, k2, k3, k4 = 0, 0, 0, 0
    result = []
    while t <  tmax + deltat:
        result.append([y, y- np.exp(-1*t)])
        k1 = lam* y
        k2 = lam* (y + deltat* k1/2)
        k3 = lam* (y + deltat*k2/2)
        k4 = lam* (y + deltat*k3)
        y += deltat * (k1 + 2* k2 + 2* k3 + k4)/6
        t += deltat
    return result

def Vn(t, Wv):
    return 10*math.sin(Wv*t)

def rrz_method(deltat, i0, q0, Wv, tmax, tmin, L, C, R):
    t = tmin
    i = i0
    q = q0
    result = []
    while t < tmax + deltat:
        result.append([i, q])
        k1_q = i
        k1_i = Vn(t, Wv)/L - q/(L*C) - R*i/L
        k2_q = i + deltat*k1_i/2
        k2_i = Vn(t+deltat/2, Wv)/L - (q + deltat*k1_q/2)/(L*C) - R*(i + deltat*k1_i/2)/L
        k3_q = i + deltat*k2_i/2
        k3_i = Vn(t+deltat/2, Wv)/L - (q + deltat*k2_q/2)/(L*C) - R*(i + deltat*k2_i/2)/L
        k4_q = i + deltat*k3_i/2
        k4_i = Vn(t+deltat, Wv)/L - (q + deltat*k3_q)/(L*C) - R*(i + deltat*k3_i)/L

        i += deltat /6*(k1_i + 2*k2_i + 2*k3_i + k4_i)
        q += deltat /6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
        t += deltat
    return result


deltat = [0.01, 0.1, 1.0]
filename = ["euler1", "euler2", "euler3"]
filename_error = ["euler1_error", "euler2_error", "euler3_error"]
x=[]
x.append([i/int(deltat[0]**-1) for i in range(0, 502)])
x.append([i/int(deltat[1]**-1) for i in range(0, 52)])
x.append([i/int(deltat[2]**-1) for i in range(0, 6)])

# Metoda Eulera
y = [list(zip(*euler_method(deltat[i])))[0] for i in range(3)]
y_err = [list(zip(*euler_method(deltat[i])))[1] for i in range(3)]
draw(x, y, 'euler_plot.png', 'z.1 - Metoda Eulera - ', y_err, 1, 'euler_error.png')

# Metoda RK2
y = [list(zip(*rk2_method(deltat[i])))[0] for i in range(3)]
y_err = [list(zip(*rk2_method(deltat[i])))[1] for i in range(3)]
draw(x, y, 'rk2_plot.png', 'z.1 - Metoda RK2 - ', y_err, 1, 'rk2_error.png', 1)

# Metoda RK4
y = [list(zip(*rk4_method(deltat[i])))[0] for i in range(3)]
y_err = [list(zip(*rk4_method(deltat[i])))[1] for i in range(3)]
draw(x, y, 'rk4_plot.png', 'z.1 - Metoda RK4 - ', y_err, 1, 'rk4_error.png', 1)

# RRZ 2 rzedu
R = 100
L = 0.1
C = 0.001
i_0 = 0
q_0 = 0
delt = math.pow(10, -4)
omega_0 = 100
omega_v = [0.5*omega_0, 0.8*omega_0, 1.0*omega_0, 1.2*omega_0]
T_0 = math.pi*2/omega_0
t_max = 4*T_0
t_min = 0

x = ([i for i in range(0, int(t_max)+2)])

x = np.linspace(t_min, t_max, int(t_max/delt)+2)
rrz1 = rrz_method(delt, i_0, q_0, omega_v[0], t_max, t_min, L, C, R)
rrz2 = rrz_method(delt, i_0, q_0, omega_v[1], t_max, t_min, L, C, R)
rrz3 = rrz_method(delt, i_0, q_0, omega_v[2], t_max, t_min, L, C, R)
rrz4 = rrz_method(delt, i_0, q_0, omega_v[3], t_max, t_min, L, C, R)
plt.figure()
plt.plot(x, list(zip(*rrz1))[0], color='blue', linewidth = 2)
plt.plot(x, list(zip(*rrz2))[0], color='red', linewidth = 2)
plt.plot(x, list(zip(*rrz3))[0], color='green', linewidth = 2)
plt.plot(x, list(zip(*rrz4))[0], color='lightblue', linewidth = 2)
plt.legend(['0.5 omega_0', '0.8 omega_0', '1.0 omega_0', '1.2 omega_0'], loc = 'upper right')
plt.ylabel('I(t)')
plt.xlabel('t')
plt.title('z.4 - Metoda RK4, I(t)')
plt.savefig('rrz_i.png')

plt.figure()
plt.plot(x, list(zip(*rrz1))[1], color='blue', linewidth = 2)
plt.plot(x, list(zip(*rrz2))[1], color='red', linewidth = 2)
plt.plot(x, list(zip(*rrz3))[1], color='green', linewidth = 2)
plt.plot(x, list(zip(*rrz4))[1], color='lightblue', linewidth = 2)
plt.legend(['0.5 omega_0', '0.8 omega_0', '1.0 omega_0', '1.2 omega_0'], loc = 'upper right')
plt.ylabel('Q(t)')
plt.xlabel('t')
plt.title('z.4 - Metoda RK4, Q(t)')
plt.savefig('rrz_q.png')