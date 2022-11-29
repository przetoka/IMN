import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
import numpy as np

sum = ["sum16.txt", "sum8.txt", "sum4.txt", "sum2.txt", "sum1.txt"]

plt.xlabel('it')
plt.ylabel('S')
plt.title('S(it)')

for i in range(5):
    data2  = np.loadtxt("data/" + sum[i],usecols=(0, 1)) 
    x = data2[:, 0]
    y = data2[:, 1]
    plt.plot(x, y)
    plt.legend(['k = 16', 'k = 8', 'k = 4', 'k = 2', 'k = 1'])
    
plt.savefig("plots/" + "sum.png")
plt.clf()




k = [16, 8, 4, 2, 1]
maps = ["map16.txt", "map8.txt", "map4.txt", "map2.txt", "map1.txt"]


for i in range(5):
    data = np.loadtxt("data/" + maps[i],usecols=(0, 1, 2)) 

    x = data[:, 0]
    y = data[:, 1]
    z = (data[:, 2])

    yi = np.arange(0.0, 25.6, y[1])
    xi = np.arange(0.0, 25.6, y[1])
    
    gx,gy = np.meshgrid(xi,yi)

    gridT12 = griddata( (x, y), z, (gx, gy), method="nearest")

    plt.title('k = ' + str(k[i]))
    plt.ylabel("y")
    plt.xlabel("x")
    pa=plt.pcolormesh(xi,yi, gridT12, shading='auto')

    cbar = plt.colorbar(pad=.01)
    plt.jet()

    plt.savefig("plots/map" + str(k[i]) + ".png")
    plt.clf()