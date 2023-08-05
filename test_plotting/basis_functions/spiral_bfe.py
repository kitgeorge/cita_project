import matplotlib.pyplot as plt
import numpy as np
import math

N_R = 100
# N_phi = 360

R_Ka = 20
R_values = np.array([i/N_R*R_Ka for i in range(N_R)])
# phi_values = np.array([i/N_phi*2*math.pi for i in range(N_phi)])

density = np.loadtxt("../../data/basis_functions/spiral_density.csv")
trunc_density = np.loadtxt("../../data/basis_functions/spiral_trunc_density.csv")
trunc_potential = np.loadtxt("../../data/basis_functions/spiral_trunc_potential.csv")
trunc_force = np.loadtxt("../../data/basis_functions/spiral_trunc_force.csv")

# density = density.reshape(N_R, N_phi)
# trunc_density = trunc_density.reshape(N_R, N_phi)
# trunc_potential = trunc_potential.reshape(N_R, N_phi)
# trunc_force = trunc_force.reshape(N_R, N_phi, 2)

plt.plot(R_values, density[:], label="Density")
plt.plot(R_values, trunc_density[:], label="Truncated Density")
plt.xlim(0, R_Ka)
plt.savefig("../../test_plots/basis_functions/spiral_trunc_density.png", dpi=500)
