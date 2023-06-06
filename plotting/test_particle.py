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
N_timesteps = 10001
data = data.reshape(N_tp, N_timesteps, 2, 2)
data = data[:, :9999, :, :]

R_0 = 8*3.086e19
R_min = 0.9*R_0
R_max = 1.1*R_0

actions_data = np.loadtxt("../data/test_particle_1_aa.csv")
actions_data = actions_data.reshape(N_tp, N_timesteps - 1, 2, 2)

kpc = 3.086e19
Myr = 1e6*365*24*60*60
Omega = 220e3/(8*kpc)
Omega_p = (1 + math.sqrt(2)/2)*Omega
integration_time = 50*2*math.pi/Omega
timestep = integration_time/1e4
time_values = np.array([i*timestep for i in range(N_timesteps - 1)])
slow_angles = np.array([actions_data[i, :, 1, 0] + 2*(actions_data[i, :, 1, 1] 
                        - Omega_p*time_values) for i in range(N_tp)]) % (2*math.pi)

apocentre_indices = [[] for i in range(N_tp)]

for i in range(N_tp):
    for j in range(1, N_timesteps - 3):
        if(data[i, j, 0, 0] > data[i, j - 1, 0, 0] 
           and data[i, j, 0, 0] > data[i, j + 1, 0, 0]):
            apocentre_indices[i].append(j)


# for i in range(N_tp):
#     plt.figure()
#     plt.plot(data[i, :, 0, 0])
#     # plt.plot(data[i, :, 0, 1]*1e20)
#     plt.plot(data[i, :, 0, 0]*data[i, :, 1, 1]/220e3)
#     plt.savefig("../plots/test_particle/tp_{}.png".format(i), dpi=500)
#     plt.close()
for j in apocentre_indices[0]:
    print(slow_angles[0, j], actions_data[i, j, 0, 1]/2)
for i in range(N_tp):
    plt.figure()
    plt.scatter([slow_angles[i, j] for j in apocentre_indices[i]], 
                [actions_data[i, j, 0, 1]/2 for j in apocentre_indices[i]],
             c=[time_values[j]/integration_time for j in apocentre_indices[i]], cmap="jet", marker=".")
    plt.xlim(0, 2*math.pi)
    plt.colorbar()
    # plt.plot(time_values, actions_data[i, :, 0, 0])
    plt.savefig("../plots/test_particle_aa/tp_{}.png".format(i), dpi=500)
    plt.close()


# for i in range(100):
#     figure = plt.figure()
#     ax = figure.add_subplot(111, polar='True')
#     c = ax.scatter(data[:, 100*i, 0, 1], data[:, 100*i, 0, 0])
#     plt.savefig("../plots/test_particle/t={}.png".format(i), dpi=500)
#     plt.close()
