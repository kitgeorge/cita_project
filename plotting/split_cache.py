import numpy as np
import pandas as pd

k_max = 10
n_max = 50
l_max = 50
N_R_tabulated = 10000
table_shape = (k_max + 1, (n_max + 1)*(l_max + 1)*N_R_tabulated)

print("Reading u values")
u_values = pd.read_csv("../cache/basis_functions/u_values.csv").to_numpy()[:, 0]
print("Reading u prime values")
u_prime_values = pd.read_csv("../cache/basis_functions/u_prime_values.csv").to_numpy()[:, 0]
print("Reading d values")
d_values = pd.read_csv("../cache/basis_functions/d_values.csv").to_numpy()[:, 0]

u_values = u_values.reshape(table_shape)
u_prime_values = u_prime_values.reshape(table_shape)
d_values = d_values.reshape(table_shape)

for i in range(k_max + 1):
    np.savetxt("../data/split_cache/u/k={}.csv".format(i), u_values[i, :])
    np.savetxt("../data/split_cache/u_prime/k={}.csv".format(i), u_prime_values[i, :])
    np.savetxt("../data/split_cache/d/k={}.csv".format(i), d_values[i, :])
