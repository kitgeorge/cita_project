import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 1.0),
    "axes.facecolor": (1.0, 1.0, 1.0, 1.0),
    "savefig.facecolor":(1.0, 1.0, 1.0, 1.0),
})
import sys
sys.path.append('/home/kit/Documents/cita_project/')
from py_src.disc_plotter import DiscPlotter
import math
import matplotlib
matplotlib.use("Agg")

data = np.loadtxt("../data/debug.csv")
N_R = 100
N_phi = 100
N_t = 100
data = data.reshape(N_R, N_phi, N_t, 2)

R_max = 2*8
R_interval = R_max/N_R
R_min = 0#50*R_interval

for k in range(1, N_t):
    print(k)
    for i in range(N_R):
        for j in range(N_phi):
            if(data[i, j, k, 0] != 0):
                print("x, {}, {}, {}".format(i, j, k))
            if(data[i, j, k, 1] != 0):
                print("x, {}, {}, {}".format(i, j, k))
print("A")
plotter = DiscPlotter(R_min, R_max, R_interval, N_phi)
for i in range(N_t):
    plotter.plotColor(data[:, :, i, 0], "../plots/debug/debug_{}.png".format(i))

# for i in range(N_t):
#     plt.figure()
#     plt.plot(data[:, 0, i, 0])
#     plt.savefig("../plots/debug/debug_{}.png".format(i))
#     plt.close()