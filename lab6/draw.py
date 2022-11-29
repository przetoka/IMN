import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from scipy.interpolate import griddata
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

maps = ["5_a.txt", "5_b.txt", "5_c.txt", "6_a.txt", "6_b.txt", "6_c.txt"]
maps_png = ["5_a.png", "5_b.png", "5_c.png", "6_a.png", "6_b.png", "6_c.png"]


for i in range(6):
    data = np.loadtxt("data/" + maps[i],usecols=(0, 1, 2)) 

    x = data[:, 0]
    y = data[:, 1]
    z = (data[:, 2])

    yi = np.arange(0.0, x[-1] + 0.1, x[1])
    xi = np.arange(0.0, x[-1] + 0.1, x[1])
    
    gx,gy = np.meshgrid(xi,yi)

    gridT12 = griddata( (x, y), z, (gx, gy), method="nearest")

    plt.title(' ')
    plt.ylabel("y")
    plt.xlabel("x")
    pa=plt.pcolormesh(xi,yi, gridT12, shading='auto', cmap = 'bwr')
    if (i >= 3): 
        plt.clim(-0.8,0.8)
    else:
        plt.clim(-10,10)
    cbar = plt.colorbar(pad=.01)

    plt.savefig("plots/map" + maps_png[i])
    plt.clf()