import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata

gl = ["global_sum1.txt", "global_sum2.txt"]
plt.xlabel('nr iteracji')
plt.ylabel('S')
plt.title('Relaksacja globalna')

for i in gl:
    fileout2 = i
    data2  = np.loadtxt("data/" + fileout2,usecols=(0, 1)) 
    x = data2[:, 0]
    y = data2[:, 1]
    plt.xlim([1, 100000])
    plt.plot(x, y)
    plt.legend(['w = 0.6', 'w = 1.0'])
    plt.xscale('log')
    
plt.savefig("plots/" + "global.png")
plt.clf()


local = ["local_sum1.txt", "local_sum2.txt", "local_sum3.txt", "local_sum4.txt"]
plt.xlabel('nr iteracji')
plt.ylabel('S')
plt.title('Relaksacja lokalna')

for i in local:
    fileout2 = i
    data2  = np.loadtxt("data/" + fileout2,usecols=(0, 1)) 
    x = data2[:, 0]
    y = data2[:, 1]
    plt.plot(x, y)
    plt.xlim([1, 100000])
    plt.legend(['w = 1.0', 'w = 1.4', 'w = 1.8', 'w = 1.9'])
    plt.xscale('log')
    
plt.savefig("plots/" + "local.png")
plt.clf()

ax=plt.subplot(111)
yi = np.arange(0.0, 10.0) # pprzesuniety zakres
xi = np.arange(0.0, 15.0)

gx,gy = np.meshgrid(xi,yi)

err = ["global_err1.txt", "global_err2.txt"]
for i in err:
    fileout2 = i
    data2  = np.loadtxt("data/" + fileout2,usecols=(0, 1, 2)) 

    x = data2[:, 0]
    y = data2[:, 1]
    z = (data2[:, 2])

    gridT12 = griddata( (x, y), z, (gx, gy), method="nearest")
    cmap=cm.summer
    
    plt.ylabel("y")
    plt.xlabel("x")
    pa=plt.pcolormesh(xi,yi, gridT12, cmap=cmap, rasterized='true', shading='gouraud')

    if i == err[0]:
        plt.title('Relaksacja globalna w = 0.6')
    else:
        plt.title('Relaksacja globalna w = 1.0')

    cbar = plt.colorbar(pad=.01)
    plt.viridis()

    plt.savefig("plots/" + i.split(".")[0] + ".png")
    plt.clf()

maps = ["global_map1.txt", "global_map2.txt"]
for i in maps:
    fileout2 = i
    data2  = np.loadtxt("data/" + fileout2,usecols=(0, 1, 2)) 

    y = data2[:, 0]
    x = data2[:, 1]
    z = (data2[:, 2])

    gridT12 = griddata( (y, x), z, (gx, gy), method="nearest")
    cmap=cm.summer
    
    plt.ylabel("y")
    plt.xlabel("x")
    pa=plt.pcolormesh(xi,yi, gridT12, cmap=cmap, rasterized='true', shading='gouraud')
    plt.viridis()
    if i == maps[0]:
        plt.title('Relaksacja globalna w = 0.6')
    else:
        plt.title('Relaksacja globalna w = 1.0')
    cbar = plt.colorbar(pad=.01)
    
    plt.savefig("plots/" + i.split(".")[0] + ".png")
    plt.clf()
