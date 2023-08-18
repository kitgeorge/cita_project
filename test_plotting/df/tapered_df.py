import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../../test_data/df/tapered_df.csv")
N_E = 1000
N_L = 1000
data = data.reshape((N_E, N_L))

plt.imshow(data, origin="lower")
plt.savefig("../../test_plots/df/tapered_df.png", dpi=500)