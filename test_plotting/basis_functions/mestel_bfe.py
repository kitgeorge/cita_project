import matplotlib.pyplot as plt
import numpy as np

trunc_density = np.loadtxt("../../test_data/basis_functions/trunc_density.csv")
trunc_potential = np.loadtxt("../../test_data/basis_functions/trunc_potential.csv")
trunc_force = np.loadtxt("../../test_data/basis_functions/trunc_force.csv")

mestel_density = np.loadtxt("../../test_data/basis_functions/mestel_density.csv")
mestel_potential = np.loadtxt("../../test_data/basis_functions/mestel_potential.csv")
mestel_force = np.loadtxt("../../test_data/basis_functions/mestel_force.csv")

N_R = 1000
R_Ka = 20
R_values = np.array([i/N_R*R_Ka for i in range(N_R)])

trunc_force.reshape(N_R, 2)
mestel_force.reshape(N_R, 2)


plt.plot(R_values, mestel_density, label="Mestel")
plt.plot(R_values, trunc_density, label="Truncated")
plt.xlim(0, 20)
plt.ylim(0, 3e9)
plt.legend()
plt.savefig("../../test_plots/basis_functions/trunc_density.png", dpi=500)
plt.close()

plt.plot(R_values, mestel_potential, label="Mestel")
plt.plot(R_values, trunc_potential, label="Truncated")
plt.xlim(0, 20)
plt.legend()
plt.savefig("../../test_plots/basis_functions/trunc_potential.png", dpi=500)
plt.close()

plt.plot(R_values, mestel_force[:, 0], label="Mestel")
plt.plot(R_values, trunc_density[:, 0], label="Truncated")
plt.legend()
plt.savefig("../../test_plots/basis_functions/trunc_f_R.png", dpi=500)
plt.close()

plt.plot(R_values, mestel_force[:, 1], label="Mestel")
plt.plot(R_values, trunc_density[:, 1], label="Truncated")
plt.legend()
plt.savefig("../../test_plots/basis_functions/trunc_f_phi.png", dpi=500)
plt.close()

