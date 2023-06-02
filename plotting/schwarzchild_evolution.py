import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 1.0),
    "axes.facecolor": (1.0, 1.0, 1.0, 1.0),
    "savefig.facecolor":(1.0, 1.0, 1.0, 1.0),
})
import math

data = np.loadtxt("../data/test_particle_1.csv")
N_tp = 1000
N_timesteps = 10001
data = data.reshape(N_tp, N_timesteps, 2, 2)

R_0 = 8*3.086e19
R_min = 0.9*R_0
R_max = 1.1*R_0


figure = plt.figure()
ax = figure.add_subplot(111, polar='True')
c = ax.scatter(data[:, 999, 0, 1], data[:, 999, 0, 0])
plt.show()
