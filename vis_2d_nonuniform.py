import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import sys
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap


cmap_name = 'coolwarm'
filename = 'data/temp/000000.out.txt'
plot_dat = 'temp'
max_lvl = 9

max_lvl   = int(sys.argv[1])
filename  = sys.argv[2]
imgfile   = sys.argv[3]
cmap_name = sys.argv[4]
plot_dat  = sys.argv[5]

cells = np.loadtxt(filename, delimiter=',')
# print(cells)

X = cells[:,1]
Y = cells[:,2]
Z = cells[:,3]
L = cells[:,0]

# define grid
xi = np.linspace(0,2**max_lvl,2**max_lvl)
yi = np.linspace(0,2**max_lvl,2**max_lvl)


cmp = cm.get_cmap(cmap_name)
tcks = []
cbar_label = ''

zi = griddata((Y, X), Z, (xi[None,:], yi[:,None]), method='linear')
if plot_dat == 'lvls':
    zi = griddata((Y, X), L, (xi[None,:], yi[:,None]), method='nearest')
    viridis = cm.get_cmap(cmap_name, 12)
    newcolors = viridis(np.linspace(0, 1, int(L.max() - L.min()) + 1))
    cmp = ListedColormap(newcolors)
    tcks = [i for i in range(int(L.min()), int(L.max()) + 1)]
    cbar_label = 'grid level'

elif plot_dat == 'procs':
    n_procs = int(sys.argv[6])

    viridis = cm.get_cmap(cmap_name, 100)
    newcolors = viridis(np.linspace(0, 1, n_procs))
    cmp = ListedColormap(newcolors)
    tcks = [i for i in range(n_procs)]
    cbar_label = 'processor number'

    C = []
    sum = 0
    for p in range(n_procs-1):
        C += [p for i in range (len(X)//n_procs)]
        sum += len(X)//n_procs
    C += [n_procs for i in range (len(X) - sum)]
    C = np.array(C)
    print(len(X), len(Y), len(C))
    zi = griddata((Y, X), C, (xi[None,:], yi[:,None]), method='nearest')
    

pl = plt.pcolor(xi, yi, zi, cmap=cmp)
cbar = plt.colorbar(pl, ticks=tcks)
if cbar_label != '':
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbar_label, rotation=270)
# plt.show()
plt.savefig(imgfile)


# plt.scatter(Y, X, c=Z, cmap=plt.get_cmap(cmap_name))
# plt.tripcolor(X, Y, L, cmap="RdBu_r")
# ax2.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
# plt.show()

# fig, ax = plt.subplots()
# p = ax.pcolor(X, Y, Z, cmap=cm.RdBu, vmin=abs(Z).min(), vmax=abs(Z).max())
# cb = fig.colorbar(p)
# plt.show()

