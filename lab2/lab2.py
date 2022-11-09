import matplotlib.pyplot as plt
import math
import numpy as np

def draw(x, y, filename, title):
    plt.plot(x, y[0], color = 'darkorange')
    plt.plot(x, y[1], color = 'royalblue')
    plt.title(title)
    plt.legend(['u(t)', 'v(t)'])
    plt.xlabel('t')
    plt.ylabel('u(t), v(t)')
    plt.savefig(filename)
    plt.clf()

def picards_method(beta, N, gamma, t_max, deltat, u_0, TOL, mu_max):
    U = [u_0]
    alpha = beta*N - gamma
    t = 0.0
    u_mu = u_0
    for i in range(int(t_max/deltat)):
        mu = 0
        while mu <= mu_max:
            temp = U[i] + deltat*((alpha*U[i] - beta*(U[i]**2)) + (alpha*u_mu - beta*(u_mu**2)))/2
            if math.fabs(temp - u_mu) < TOL:
                break
            u_mu = temp
            mu +=1
        U.append(u_mu)
        t +=deltat
    return U

def newtons_iteration(beta, N, gamma, t_max, deltat, u_0, TOL, mu_max):
    U = [u_0]
    alpha = beta*N - gamma
    t = 0.0
    u_mu = u_0
    for i in range(int(t_max/deltat)):
        mu = 0
        while mu <= mu_max:
            temp = u_mu - (u_mu - U[i] - deltat*((alpha*U[i] - beta*(U[i]**2)) + (alpha*u_mu - beta*(u_mu**2)))/2)/(1-deltat*(alpha - 2*beta*u_mu)/2)
            if math.fabs(temp - u_mu) < TOL:
                break
            u_mu = temp
            mu +=1
        U.append(u_mu)
        t +=deltat
    return U

def f(t, u):
    return (beta*N - gamma)*u - beta*u*u

def rk2_method(beta, N, gamma, t_max, deltat, u_0, TOL, mu_max, a, b, c):
    U = [u_0]
    U1 = U2 = U[0]
    alpha = beta*N - gamma
    t = 0.0
    m = [[0, 0], [0, 0]]
    for i in range(int(t_max/deltat)):
        mu = 0
        U1 = U[i]
        U2 = U[i]

        F1 = U1 - U[i] - deltat*(a[0][0]*(alpha*U1 - beta*U1**2) + a[0][1]*(alpha*U2 - beta*U2**2))
        F2 = U2 - U[i] - deltat*(a[1][0]*(alpha*U1 - beta*U1**2) + a[1][1]*(alpha*U2 - beta*U2**2))

        m[0][0] = 1 - deltat*a[0][0]*(alpha - 2*beta*U1)
        m[0][1] =  - deltat*a[0][1]*(alpha - 2*beta*U2)
        m[1][0] =  - deltat*a[1][0]*(alpha - 2*beta*U1)
        m[1][1] = 1 - deltat*a[1][1]*(alpha - 2*beta*U2)

        deltaU1 = (F2 * m[0][1] - F1 * m[1][1])/(m[0][0] * m[1][1] - m[0][1] * m[1][0])
        deltaU2 = (F1 * m[1][0] - F2 * m[0][0])/(m[0][0] * m[1][1] - m[0][1] * m[1][0])

        u_mu = U[i]
        while mu <= mu_max:
            mu +=1
            temp = U[i] + deltat*(f(t+(1/2 - math.sqrt(3)/6)*deltat, U1)/2 + f(t+(1/2 + math.sqrt(3)/6)*deltat, U2)/2)
            U1 += deltaU1
            U2 += deltaU2
            if math.fabs(temp - u_mu) < TOL:
                break
            u_mu = temp
        U.append(u_mu)
        t +=deltat
    return U

beta = 0.001
N = 500
gamma = .1
t_max = 100
deltat = .1
u_0 = 1
TOL = 10**(-6)
mu_max = 20
x = np.linspace(0.0, t_max, int(t_max/deltat)+1)
a = [[1/4, 1/4 - math.sqrt(3)/6], [1/4 + math.sqrt(3)/6, 1/4]]
b = [1/2, 1/2]
c = [1/2 - math.sqrt(3)/6, 1/2 + math.sqrt(3)/6]

# Metoda Picarda
picard = picards_method(beta, N, gamma, t_max, deltat, u_0, TOL, mu_max)
y = [picard, [N - picard[i] for i in range(0, int(t_max/deltat)+1)]]
draw(x, y, 'picard.png', 'Metoda Picarda')

# Metoda Newtona
newton = newtons_iteration(beta, N, gamma, t_max, deltat, u_0, TOL, mu_max)
y = [newton, [N - newton[i] for i in range(0, int(t_max/deltat)+1)]]
draw(x, y, 'newton.png', 'Iteracja Newtona')

# Metoda RK2
rk2 = rk2_method(beta, N, gamma, t_max, deltat, u_0, TOL, mu_max, a, b, c)
y = [rk2, [N - rk2[i] for i in range(0, int(t_max/deltat)+1)]]
draw(x, y, 'rk2.png', 'Niejawna RK2')