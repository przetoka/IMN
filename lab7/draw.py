import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from scipy.interpolate import griddata
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

maps = ["psi_zeta_-1000.txt", "psi_zeta_-4000.txt", "psi_zeta_4000.txt"]
maps_png = ["psi_-1000.png", "zeta_-1000.png", "psi_-4000.png", "zeta_-4000.png", "psi_4000.png", "zeta_4000.png"]
titles = ["Q = -1000, ψ(x, y)", "Q = -1000, ζ(x, y)", "Q = -4000, ψ(x, y)", "Q = -4000, ζ(x, y)", "Q = 4000, ψ(x, y)", "Q = 4000, ζ(x, y)"]
c = [-55, -50, -218, -202, 202, 218]

for i in range(3):
    data = np.loadtxt("data/" + maps[i],usecols=(0, 1, 2, 3)) 

    x = data[:, 0]
    y = data[:, 1]
    psi = (data[:, 2])
    zeta = (data[:, 3])


    levels = np.linspace(psi.min(), psi.max(), 600)

    plt.title(titles[i])
    plt.ylabel("y")
    plt.xlabel("x")
    plt.tricontour(x, y, psi, levels)
    plt.set_cmap('gnuplot')
    cbar = plt.colorbar(pad=.01)
    plt.clim(c[2*i], c[2*i + 1])
    plt.savefig("plots/" + maps_png[2*i])
    plt.clf()

    levels = np.linspace(zeta.min(), zeta.max(), 200)
    
    plt.title(titles[i + 1])
    plt.ylabel("y")
    plt.xlabel("x")
    plt.tricontour(x, y, zeta, levels)
    plt.set_cmap('gnuplot')
    cbar = plt.colorbar(pad=.01)

    plt.savefig("plots/" + maps_png[2*i + 1])
    plt.clf()


maps = ["u_v_-1000.txt", "u_v_-4000.txt"]
maps_png = ["u_-1000.png", "v_-1000.png", "u_-4000.png", "v_-4000.png"]
titles = ["Q = -1000, u(x, y)", "Q = -1000, v(x, y)", "Q = -4000, u(x, y)", "Q = -4000, v(x, y)"]


for i in range(len(maps)):
    data = np.loadtxt("data/" + maps[i],usecols=(0, 1, 2, 3)) 

    x = data[:, 0]
    y = data[:, 1]
    u = (data[:, 2])
    v = (data[:, 3])

    yi = np.arange(0.0, y[-1], y[1])
    xi = np.arange(0.0, x[-1], y[1])
    
    gx,gy = np.meshgrid(xi,yi)

    U = griddata( (x, y), u, (gx, gy), method="nearest")
    V = griddata( (x, y), v, (gx, gy), method="nearest")

    plt.title(titles[i])
    plt.ylabel("y")
    plt.xlabel("x")
    pa=plt.pcolormesh(xi,yi, U, shading='auto')
    plt.jet()
    cbar = plt.colorbar(pad=.01)

    plt.savefig("plots/" + maps_png[2*i])
    plt.clf()

    plt.title(titles[i + 1])
    plt.ylabel("y")
    plt.xlabel("x")
    pa=plt.pcolormesh(xi,yi, V, shading='auto')
    plt.jet()
    cbar = plt.colorbar(pad=.01)

    plt.savefig("plots/" + maps_png[2*i + 1])
    plt.clf()