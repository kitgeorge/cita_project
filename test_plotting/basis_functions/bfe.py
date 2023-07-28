import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../../test_data/basis_functions/basis_functions.csv")
N_n = 6
N_k = 6
N_R = 100
data = data.reshape(N_n, N_k, 2, N_R)
R_values = np.array([i/N_R for i in range(N_R)])


fig, axs = plt.subplots(4, 3, figsize=(10, 10))

for i in range(N_n):
    for j in range(N_k):
        if(i < 3):
            r = 1
            c = i
        else:
            r = 3
            c = i - 3
        axs[r, c].plot(R_values, data[i, j, 0, :], label="{}".format(j))
    for j in range(N_k):
        if(i < 3):
            r = 0
            c = i
        else:
            r = 2
            c = i - 3
        axs[r, c].plot(R_values, data[i, j, 1, :], label="{}".format(j))
plt.legend()

plt.savefig("../../test_plots/basis_functions/basis_functions_combined.png", dpi=500)