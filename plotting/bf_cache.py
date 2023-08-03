import numpy as np
import matplotlib.pyplot as plt

k_max = 10
n_max = 50
l_max = 50
N_R_tabulated = 10000
table_shape = (k_max + 1, n_max + 1, l_max + 1, N_R_tabulated)

u_values = np.loadtxt("../cache/basis_functions/u_values.csv")
u_prime_values = np.loadtxt("../cache/basis_functions/u_prime_values.csv")
d_values = np.loadtxt("../cache/basis_functions/d_values.csv")
u_values = u_values.reshape(table_shape)
u_prime_values = u_prime_values.reshape(table_shape)
d_values = d_values.reshape(table_shape)

k = 4
l = 2

for i in range(n_max + 1):
    plt.figure()
    plt.plot(u_values[k, l, i, :])
    plt.savefig("../plots/bf_cache/u/k=4l=2/n={}.png".format(i), dpi=500)
    plt.close()

    plt.plot(u_prime_values[k, l, i, :])
    plt.savefig("../plots/bf_cache/u_prime/k=4l=2/n={}.png".format(i), dpi=500)
    plt.close()

    plt.plot(d_values[k, l, i, :])
    plt.savefig("../plots/bf_cache/d/k=4l=2/n={}.png".format(i), dpi=500)
    plt.close()
