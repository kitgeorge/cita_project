import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 1.0),
    "axes.facecolor": (1.0, 1.0, 1.0, 1.0),
    "savefig.facecolor":(1.0, 1.0, 1.0, 1.0),
})
import math

data = np.loadtxt("../data/test_particle_1.csv")
N_tp = 36
N_timesteps = 100001
data = data.reshape(N_tp, N_timesteps, 2, 2)

R_0 = 8*3.086e19
R_min = 0.9*R_0
R_max = 1.1*R_0



for i in range(N_tp):
    plt.figure()
    plt.plot(data[i, :, 0, 0])
    # plt.plot(data[i, :, 0, 1]*1e20)
    plt.plot(data[i, :, 0, 0]*data[i, :, 1, 1]/220e3)
    plt.savefig("../plots/test_particle/tp_{}.png".format(i), dpi=500)
    plt.close()

# for i in range(100):
#     figure = plt.figure()
#     ax = figure.add_subplot(111, polar='True')
#     c = ax.scatter(data[:, 100*i, 0, 1], data[:, 100*i, 0, 0])
#     plt.savefig("../plots/test_particle/t={}.png".format(i), dpi=500)
#     plt.close()
