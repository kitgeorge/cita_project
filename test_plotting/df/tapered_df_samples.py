import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/kit/Documents/cita_project/')
from py_src.disc_plotter import DiscPlotter

samples = np.loadtxt("../../test_data/df/tapered_df_sample.csv")
N_samples = 10000
samples = samples.reshape(N_samples, 2)

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.scatter(samples[:, 1], samples[:, 0], marker=".", s=0.001)
plt.savefig("../../test_plots/df/tapered_df_samples.png", dpi=5000)