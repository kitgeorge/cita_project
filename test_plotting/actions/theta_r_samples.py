import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../../test_data/actions/theta_r_coords.csv")

N_angles = 100
data = data.reshape(N_angles, 2)

colors = [i/N_angles for i in range(N_angles)]

plt.scatter(data[:, 0], data[:, 1], c=colors, cmap="turbo")
plt.xlabel(r'$R/\mathrm{kpc}$')
plt.ylabel(r'$v_R/\mathrm{kms}^{-1}$')

plt.savefig("../../test_plots/actions/theta_r_samples.png", dpi=500)