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

data = np.loadtxt("../data/test_particle_dehnen_aa.csv")
N_tp = 1000
N_timesteps = 500
data = data.reshape(N_tp, N_timesteps, 2, 2)

raw_data = np.loadtxt("../data/test_particle_dehnen.csv")
raw_data = raw_data.reshape(N_tp, N_timesteps, 2, 2)


for i in range(N_timesteps):
    print(i)
    plt.figure()
    plt.scatter(data[:, i, 0, 1], data[:, i, 0, 0], marker=".")
    plt.savefig("../plots/weak_transients/exp_action/timestep={}.png".format(i), dpi=500)
    plt.close("all")

for i in range(N_timesteps):
    print(i)
    plt.figure()
    plt.hist(data[:, i, 0, 0])
    plt.savefig("../plots/weak_transients/exp_action_hists/timestep={}.png".format(i), dpi=500)
    plt.close("all")

for i in range(N_timesteps):
    print(i)
    plt.figure()
    plt.scatter(data[:, i, 0, 1], data[:, i, 0, 0] - data[:, i - 1, 0, 0], marker=".")
    plt.savefig("../plots/weak_transients/exp_action_diff/timestep={}.png".format(i), dpi=500)
    plt.close("all")

for i in range(N_timesteps):
    print(i)
    plt.figure()
    plt.hist(data[:, i, 0, 0] - data[:, 0, 0, 0])
    plt.savefig("../plots/weak_transients/exp_action_diff_hists/timestep={}.png".format(i), dpi=500)
    plt.close("all")

for i in range(N_tp):
    print(i)
    plt.figure()
    plt.plot(data[i, :, 0, 1], data[i, :, 0, 0], alpha=0.3)
    plt.scatter(data[i, :, 0, 1], data[i, :, 0, 0], c=[i/N_timesteps for i in range(N_timesteps)], cmap="jet", marker=".")
    plt.savefig("../plots/weak_transients/exp_action_trajectories/tp_{}.png".format(i), dpi=500)
    plt.close("all")

# for i in range(N_timesteps):
#     print(i)
#     plt.figure()
#     plt.scatter(raw_data[:, i, 0, 1], raw_data[:, i, 0, 0], marker=".")
#     plt.savefig("../plots/weak_transients/exp_action/timestep={}.png".format(i), dpi=500)
#     plt.close("all")

# for i in range(N_timesteps):
#     print(i)
#     plt.figure()
#     plt.hist(raw_data[:, i, 0, 0])
#     plt.savefig("../plots/weak_transients/exp_action_hists/timestep={}.png".format(i), dpi=500)
#     plt.close("all")
    
# for i in range(N_timesteps):
#     print(i)
#     plt.figure()
#     plt.scatter(raw_data[:, i, 0, 1], raw_data[:, i, 0, 0] - raw_data[:, i - 1, 0, 0], marker=".")
#     plt.savefig("../plots/weak_transients/exp_action_diff/timestep={}.png".format(i), dpi=500)
#     plt.close("all")

# for i in range(N_timesteps):
#     print(i)
#     plt.figure()
#     plt.hist(raw_data[:, i, 0, 0] - raw_data[:, 0, 0, 0])
#     plt.savefig("../plots/weak_transients/exp_action_diff_hists/timestep={}.png".format(i), dpi=500)
#     plt.close("all")

# for i in range(N_tp):
#     print(i)
#     plt.figure()
#     plt.plot(raw_data[i, :, 0, 1], raw_data[i, :, 0, 0], alpha=0.3)
#     plt.scatter(raw_data[i, :, 0, 1], raw_data[i, :, 0, 0], c=[i/N_timesteps for i in range(N_timesteps)], cmap="jet", marker=".")
#     plt.savefig("../plots/weak_transients/exp_action_trajectories/tp_{}.png".format(i), dpi=500)
#     plt.close("all")
