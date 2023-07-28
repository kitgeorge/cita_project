import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../../test_data/actions/separatrices.csv")
N_separatrix_points = 100
data = data.reshape(N_separatrix_points, 3)

trajectory = np.loadtxt("../../test_data/actions/libration_trajectory.csv")
N_timesteps = 1000
trajectory = trajectory.reshape(N_timesteps, 3)

plt.plot(data[:, 0], data[:, 1])
plt.plot(data[:, 0], data[:, 2])
plt.scatter(trajectory[:, 1], trajectory[:, 2], c=[i/N_timesteps for i in range(N_timesteps)])
plt.savefig("../../test_plots/actions/separatrices.png", dpi=500)

