import numpy as np
import math
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/kit/Documents/cita_project/')

sample_data = np.loadtxt("../../test_data/df/get_df_sample_e_l.csv")
sample_data = sample_data.reshape(1000, 2)

df_data = np.loadtxt("../../test_data/df/dehnen_df_test.csv")
df_data = df_data.reshape(1000, 1000)
print(np.min(df_data), np.max(df_data))


plt.scatter(sample_data[:, 0], sample_data[:, 1], marker=".")
plt.savefig("../../test_plots/df/dehnen_sample.png", dpi=500)

plt.figure()
plt.imshow(df_data, origin="lower")
plt.savefig("../../test_plots/df/dehnen_df.png", dpi=500)
