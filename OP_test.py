import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

from sys import path
path.append('/Users/tylerklimas/Desktop/Orbital Mechanics')


cb = pd.earth

if __name__ == '__main__':
    # intial conditions
    r_mag = cb['radius'] + 1500
    v_mag = np.sqrt(cb['mu']/r_mag)
    print(v_mag)

    r0 = [r_mag, 0, 0]
    v0 = [0, v_mag, 2]

    tspan = 100 * 60

    dt = 100
    op = OP(r0, v0, tspan, dt, cb)
    op.propogate_orbit()
    op.plot_3d(show_plot=True)
