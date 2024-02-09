import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import OrbitTools as t

from math import sqrt

from sys import path
path.append('/Users/tylerklimas/Desktop/Orbital Mechanics')

tspan = 3600 * 24
dt = 100

if __name__ == '__main__':
    r_mag = pd.earth['radius'] + 400
    v_mag = sqrt(pd.earth['mu']/r_mag)

    r0 = np.array([r_mag, 0, 0])
    v0 = np.array([0, v_mag, 0])

    state0 = np.array((r0.tolist() + v0.tolist()))

    r_mag1 = pd.earth['radius'] + 1000
    v_mag1 = sqrt(pd.earth['mu']/r_mag1) * 1.3

    r01 = np.array([r_mag1, 0, 0])
    v01 = np.array([0, v_mag1, 0.3])

    op0 = OP(state0, tspan, dt)
    # op1 = OP(r01, v01, tspan, dt)

    op0.propogate_orbit()
    # op1.propogate_orbit()

    t.plot_n_orbits([op0.rs], labels=['0', '1'], show_plot=True)
