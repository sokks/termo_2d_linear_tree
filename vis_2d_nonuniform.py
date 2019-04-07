import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

cmap_name = 'coolwarm'

cells = np.loadtxt('data/temp/000000.out.txt', delimiter=',')
print(cells)

X = cells[:,1]
Y = cells[:,2]
Z = cells[:,3]
C = np.array(['b' for i in range(len(X)//2)] + ['r' for i in range(len(X)//2)])

plt.scatter(X, Y, c=Z, cmap=plt.get_cmap(cmap_name))
plt.show()

# fig, ax = plt.subplots()
# p = ax.pcolor(X, Y, Z, cmap=cm.RdBu, vmin=abs(Z).min(), vmax=abs(Z).max())
# cb = fig.colorbar(p)
# plt.show()


alpha = 0.7
phi_ext = 2 * np.pi * 0.5

def flux_qubit_potential(phi_m, phi_p):
    return 2 + alpha - 2 * np.cos(phi_p)*np.cos(phi_m) - alpha * np.cos(phi_ext - 2*phi_p)

phi_m = np.linspace(0, 2*np.pi, 100)
phi_p = np.linspace(0, 2*np.pi, 100)
X,Y = np.meshgrid(phi_p, phi_m)
print(X)
print(Y)

Z = flux_qubit_potential(X, Y).T

fig, ax = plt.subplots()
p = ax.pcolor(X/(2*np.pi), Y/(2*np.pi), Z, cmap=cm.RdBu, vmin=abs(Z).min(), vmax=abs(Z).max())
cb = fig.colorbar(p)
plt.show()
