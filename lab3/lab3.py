import matplotlib.pyplot as plt
import math

def draw(x1, y1, x2, y2, filename, title, xlabel = 'x', ylabel = 'y', TOL1 = 2, TOL2 = 5):
    plt.plot(x1, y1, color = 'darkorange', linewidth = 2)
    plt.plot(x2, y2, color = 'blue', linewidth = 1, linestyle = 'dashed')
    plt.title(title)
    plt.legend(["TOL = 10$^{-" + str(TOL1)+ "}$", "TOL = 10$^{-" +str(TOL2) + "}$"])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(filename)
    plt.clf()

def E_x(x_1, x_2):
    return (x_2 - x_1)/(2**2 - 1)

def E_v(v_1, v_2):
    return (v_2 - v_1)/(2**2 - 1)

def rk2_method(deltat, x_n, v_n, alpha):
    k1_x = v_n
    k1_v = alpha* v_n* (1 - x_n**2) - x_n
    k2_x = v_n + deltat * k1_v
    k2_v = alpha * (1 - (x_n + deltat * k1_x)**2) * (v_n + deltat * k1_v) - (x_n + deltat * k1_x)
        
    x_new = x_n + deltat * (k1_x + k2_x)/2
    v_new = v_n + deltat * (k1_v + k2_v)/2
    return [x_new, v_new]

def f(v):
    return v

def g(x, v, alpha = 5):
    return alpha*(1 - x**2) * v - x

def F(x_k, x_n, v_k, v_n, deltat):
    return x_k - x_n - deltat*(f(v_n) + f(v_k))/2

def G(x_k, x_n, v_k, v_n, deltat):
    return v_k - v_n - deltat*(g(x_n, v_n) + g(x_k, v_k))/2

def trapezoid_method(deltat, x_n, v_n, alpha):
    x_k = x_n
    v_k = v_n
    
    a = [[1, -deltat/2], [-deltat*(-2*alpha*x_k*v_k - 1)/2, 1 - deltat*alpha*(1 - x_k**2)/2]]
    deltax = (-F(x_k, x_n, v_k, v_n, deltat) * a[1][1] + G(x_k, x_n, v_k, v_n, deltat) * a[0][1])/(a[0][0] * a[1][1] - a[0][1] * a[1][0])
    deltay = (-G(x_k, x_n, v_k, v_n, deltat) * a[0][0] + F(x_k, x_n, v_k, v_n, deltat) * a[1][0])/(a[0][0] * a[1][1] - a[0][1] * a[1][0])

    while math.fabs(deltax) > 10**(-10) and math.fabs(deltay) > 10**(-10):
        a = [[1, -deltat/2], [-deltat*(-2*alpha*x_k*v_k - 1)/2, 1 - deltat*alpha*(1 - x_k**2)/2]]
        F_temp = F(x_k, x_n, v_k, v_n, deltat)
        G_temp = G(x_k, x_n, v_k, v_n, deltat)

        deltax = (-F_temp * a[1][1] - (-G_temp) * a[0][1])/(a[0][0] * a[1][1] - a[0][1] * a[1][0])
        deltav = (-G_temp * a[0][0] - (-F_temp) * a[1][0])/(a[0][0] * a[1][1] - a[0][1] * a[1][0])

        x_k += deltax
        v_k += deltav

    return [x_k, v_k]

def time_control_algorythm(pow, method ,alpha = 5, x_0 = .01, v_0 = 0, deltat = 1, t_0 = 0, t_max = 40, S = 0.75):
    TOL = 10**(-pow)
    # r_1 = method(deltat, x_0, v_0, alpha)
    # r_2 = method(deltat*2, x_0, v_0, alpha)

    x = [x_0]
    v = [v_0]
    delta = [deltat]
    t = [t_0]

    while t[-1] <= t_max:
        r_2 = method(2*delta[-1], x[-1], v[-1], alpha)
        r_1 = method(delta[-1], x[-1], v[-1], alpha)
        r_1 = method(delta[-1], r_1[0], r_1[1], alpha)
        
        Ex = E_x(r_1[0], r_2[0])
        Ev = E_v(r_1[1], r_2[1])

        if max(math.fabs(Ex), math.fabs(Ev)) < TOL:
            x.append(r_1[0])
            v.append(r_1[1])
            delta.append(delta[-1])
            t.append(t[-1]+ delta[-1]*2)

        delta[-1] *= ((S * TOL)/max(math.fabs(Ex), math.fabs(Ev))) ** (1./3.)
    return [x, v, t, delta]


pow = [2, 5]
rk2_1 = time_control_algorythm(pow[0], rk2_method)
rk2_2 = time_control_algorythm(pow[1], rk2_method)

draw(rk2_1[2], rk2_1[0], rk2_2[2], rk2_2[0], "rk2_method_x.png", "Metoda RK2", "t", "x(t)", pow[0], pow[1])
draw(rk2_1[2], rk2_1[1], rk2_2[2], rk2_2[1], "rk2_method_v.png", "Metoda RK2", "t", "v(t)", pow[0], pow[1])
draw(rk2_1[2], rk2_1[3], rk2_2[2], rk2_2[3], "rk2_method_deltat.png", "Metoda RK2", "t", "dt(t)", pow[0], pow[1])
draw(rk2_1[0], rk2_1[1], rk2_2[0], rk2_2[1], "rk2_method_v(x).png", "Metoda RK2", "t", "v(x)", pow[0], pow[1])

trapezoid_1 = time_control_algorythm(pow[0], trapezoid_method)
trapezoid_2 = time_control_algorythm(pow[1], trapezoid_method)
draw(trapezoid_1[2], trapezoid_1[0], trapezoid_2[2], trapezoid_2[0], "trapezoid_method_x.png", "Metoda Trapezów", "t", "x(t)", pow[0], pow[1])
draw(trapezoid_1[2], trapezoid_1[1], trapezoid_2[2], trapezoid_2[1], "trapezoid_method_v.png", "Metoda Trapezów", "t", "v(t)", pow[0], pow[1])
draw(trapezoid_1[2], trapezoid_1[3], trapezoid_2[2], trapezoid_2[3], "trapezoid_method_deltat.png", "Metoda Trapezów", "t", "dt(t)", pow[0], pow[1])
draw(trapezoid_1[0], trapezoid_1[1], trapezoid_2[0], trapezoid_2[1], "trapezoid_method_v(x).png", "Metoda Trapezów", "t", "v(x)", pow[0], pow[1])
