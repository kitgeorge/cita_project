import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 1.0),
    "axes.facecolor": (1.0, 1.0, 1.0, 1.0),
    "savefig.facecolor":(1.0, 1.0, 1.0, 1.0),
})
import math
import matplotlib
matplotlib.use("Agg")

data = np.loadtxt("../data/resonance_evolution.csv")
N_J_phi = 150
N_J_R = 50
N_timesteps = 10
data = data.reshape(N_timesteps, N_J_phi, N_J_R)

timestep = 50

for i in range(N_timesteps):
    plt.imshow(np.transpose(data[i, :, :20]), origin="lower")
    plt.savefig("../plots/resonance_evolution/t={}.png".format(i*timestep),
                dpi=500)