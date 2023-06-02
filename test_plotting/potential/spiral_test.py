import numpy as np
import math
import sys
sys.path.append('/home/kit/Documents/cita_project/')
from py_src.disc_plotter import DiscPlotter


pot_data = np.loadtxt("../../test_data/potential/spiral_potential.csv")
force_data = np.loadtxt("../../test_data/potential/spiral_force.csv")
cart_force_data = np.loadtxt("../../test_data/potential/spiral_force_cart.csv")

R_max = 10
N_R_values = 100
R_interval = R_max/N_R_values
N_phi_values = 360
N_t_values = 10

pot_data = pot_data.reshape((N_R_values, N_phi_values, N_t_values))
force_data = force_data.reshape((N_R_values, N_phi_values, N_t_values, 2))
cart_force_data = cart_force_data.reshape((N_R_values, N_phi_values, N_t_values, 2))

N_arrows_R = 10
N_arrows_phi = 30

phi_values = np.array([i/N_phi_values*2*math.pi for i in range(N_phi_values)])

alt_force_data = np.zeros_like(cart_force_data)
for i in range(N_phi_values):
    alt_force_data[:, i, :, 0] = cart_force_data[:, i, :, 0]*math.cos(phi_values[i]) + cart_force_data[:, i, :, 1]*math.sin(phi_values[i])
    alt_force_data[:, i, :, 1] = cart_force_data[:, i, :, 1]*math.cos(phi_values[i]) - cart_force_data[:, i, :, 0]*math.sin(phi_values[i])

# for i in range(N_R_values):
#     for j in range(N_phi_values):
#         for k in range(N_t_values):
#             # alt_force_data[i, j, k, 0] = cart_force_data[i, j, k, 0]*math.cos(phi_values[j]) + cart_force_data[i, j, k, 1]*math.sin(phi_values[j])
#             # alt_force_data[i, j, k, 1] = cart_force_data[i, j, k, 1]*math.cos(phi_values[j]) - cart_force_data[i, j, k, 0]*math.sin(phi_values[j])
#             for l in range(2):
#                 if((alt_force_data[i, j, k, l]/force_data[i, j, k, l] - 1)**2 > 1e-4):
#                     print(l, alt_force_data[i, j, k, l], force_data[i, j, k, l])


dp = DiscPlotter(0, R_max, R_interval, N_phi_values)

for i in range(N_t_values):
    dp.plotArrowsColor(pot_data[:, :, i], alt_force_data[:, :, i, :], 
                       N_arrows_R, N_arrows_phi, 
                       "../../test_plots/potential/spiral_test/t={}.png".format(i/N_t_values))
