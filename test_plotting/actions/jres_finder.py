import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../../test_data/actions/jres_finder.csv")
N = 101
data = data.reshape(N, 2)
plt.plot(data[:, 0], data[:, 1])
plt.xlim(0, 2)
plt.ylim(0, 0.5)
plt.savefig("../../test_plots/actions/jres_finder.png", dpi=500)