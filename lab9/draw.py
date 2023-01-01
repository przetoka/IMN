#!/usr/bin/python3.9

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from scipy.interpolate import griddata
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import imageio
a = 50
maps = ["v_x_y.txt"]
maps_png = ["v_x.png", 'v_y.png']
titles = ["v_x", "v_y"]

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

data1 = np.loadtxt("data/c_1.txt", usecols=(0, 1)) 
data2 = np.loadtxt("data/c_2.txt", usecols=(0, 1)) 
data3 = np.loadtxt("data/x_sr_1.txt", usecols=(0, 1))
data4 = np.loadtxt("data/x_sr_2.txt", usecols=(0, 1))

plt.plot(data1[:, 0], data1[:, 1])
plt.plot(data2[:, 0], data2[:, 1])
plt.title("c(t_n)")
plt.ylabel("y")
plt.xlabel("x")
plt.legend(["D = 0.0", "D = 0.1"])
plt.savefig("plots/c.png")
plt.clf()

plt.plot(data3[:, 0], data3[:, 1])
plt.plot(data4[:, 0], data4[:, 1])
plt.title("x_sr(t_n)")
plt.ylabel("y")
plt.xlabel("x")
plt.legend(["D = 0.0", "D = 0.1"])
plt.savefig("plots/x_sr.png")
plt.clf()

images_0 = []
for i in range(a):
    data = np.loadtxt("data/maps/map" + str(i) + "_00.txt" ,usecols=(0, 1, 2)) 

    x = data[:, 0]
    y = data[:, 1]
    u = (data[:, 2])

    yi = np.arange(0.0, y[-1], y[1])
    xi = np.arange(0.0, x[-1], y[1])
    
    gx,gy = np.meshgrid(xi,yi)

    U = griddata( (x, y), u, (gx, gy), method="nearest")

    pa=plt.pcolormesh(xi,yi, U, shading='auto')
    plt.jet()
    cbar = plt.colorbar(pad=.01)

    plt.clim(0, 20)
    plt.savefig("plots/maps/map" + str(i) + "_00.png")
    plt.clf()
    images_0.append(imageio.v2.imread("plots/maps/map" + str(i) + "_00.png"))

images_1 = []

for i in range(a):
    data = np.loadtxt("data/maps/map" + str(i) + "_01.txt" ,usecols=(0, 1, 2)) 

    x = data[:, 0]
    y = data[:, 1]
    u = (data[:, 2])

    yi = np.arange(0.0, y[-1], y[1])
    xi = np.arange(0.0, x[-1], y[1])
    
    gx,gy = np.meshgrid(xi,yi)

    U = griddata( (x, y), u, (gx, gy), method="nearest")

    pa=plt.pcolormesh(xi,yi, U, rasterized='true', shading='gouraud')
    plt.jet()
    cbar = plt.colorbar(pad=.01)
    

    plt.clim(0, max(u))
    plt.savefig("plots/maps/map" + str(i) + "_01.png")
    plt.clf()
    images_1.append(imageio.v2.imread("plots/maps/map" + str(i) + "_01.png"))

imageio.mimsave('D_0.gif', images_0)
imageio.mimsave('D_1.gif', images_1)