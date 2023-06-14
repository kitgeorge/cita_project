import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 1.0),
    "axes.facecolor": (1.0, 1.0, 1.0, 1.0),
    "savefig.facecolor":(1.0, 1.0, 1.0, 1.0),
})
import math

class DiscPlotter:

    def __init__(self, R_min, R_max, R_interval, N_phi_intervals):
        self.N_R_intervals = int((R_max - R_min)/R_interval)
        self.N_phi_intervals = N_phi_intervals
        phi_interval = 2*math.pi/N_phi_intervals
        self.R_values = np.flip(np.array([R_max - i*R_interval for i in range(self.N_R_intervals)]))
        self.phi_values = np.array([i*phi_interval for i in range(N_phi_intervals)])
        self.X, self.Y = np.meshgrid(self.phi_values, self.R_values)

    def plotColor(self, data, path):
        # data must be indexed as (phi, R)
        figure = plt.figure()
        ax = figure.add_subplot(111, polar='True')
        pc = ax.pcolormesh(self.X, self.Y, data[:, :], shading='auto')
        figure.colorbar(pc)
        plt.savefig(path, dpi=500)
        plt.close()

    def plotArrows(self, data, N_arrows_R, N_arrows_phi, path):
        # data must match X and Y, with and extra dimension of 2
        # corresponding to f_R and f_phi
        arrow_interval_R = int(self.N_R_intervals/N_arrows_R)
        arrow_interval_phi = int(self.N_phi_intervals/N_arrows_phi)
        
        figure = plt.figure()
        ax = figure.add_subplot(111, polar='True')
        ax.quiver(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi],
                  self.Y[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi],
                  data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 0]
                  *np.cos(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi])
                  - data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 1]
                  *np.sin(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi]),
                  data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 0]
                  *np.sin(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi])
                  + data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 1]
                  *np.cos(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi]))
        plt.savefig(path, dpi=500)
        plt.close()

    def plotArrowsColor(self, data, arrows_data, N_arrows_R, N_arrows_phi, path):
        arrow_interval_R = int(self.N_R_intervals/N_arrows_R)
        arrow_interval_phi = int(self.N_phi_intervals/N_arrows_phi)
        
        figure = plt.figure()
        ax = figure.add_subplot(111, polar='True')
        ax.pcolormesh(self.X, self.Y, data[:, :], shading='auto')
        ax.quiver(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi],
                  self.Y[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi],
                  arrows_data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 0]
                  *np.cos(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi])
                  - arrows_data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 1]
                  *np.sin(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi]),
                  arrows_data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 0]
                  *np.sin(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi])
                  + arrows_data[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi, 1]
                  *np.cos(self.X[arrow_interval_R::arrow_interval_R, ::arrow_interval_phi]))
        plt.savefig(path, dpi=500)
        plt.close()
    