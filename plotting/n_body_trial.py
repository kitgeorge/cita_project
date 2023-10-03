import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/kit/Documents/cita_project/')
from py_src.disc_plotter import DiscPlotter
import math

init_theta_R_values = np.loadtxt("../data/n_body/init_theta_R_values.csv")
plt.hist(init_theta_R_values, bins=1000)
plt.savefig("../plots/init_theta_R_values.png")

coefficient_norms = np.loadtxt("../data/n_body/test_bfe_coefficient_norms.csv")
plt.plot(coefficient_norms)
plt.savefig("../plots/n_body_trial_norms.png")
# exit()

data = np.loadtxt("../data/n_body/test_trajectories_nsg.csv")
N_particles = 10000
N_timesteps = 1001
data = data.reshape(N_timesteps, N_particles, 2, 2)

for i in range(N_timesteps):
    fig = plt.figure()
    ax = fig.add_subplot(projection='polar')
    c = ax.scatter(data[i, :, 0, 1], data[i, :, 0, 0], s=0.1)
    ax.set_ylim([0, 10])
    plt.savefig("../plots/n_body_trial_nsg/t={:04d}Myr".format(i*10))
    plt.close()

