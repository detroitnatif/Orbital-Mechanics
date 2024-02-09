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

cb = pd.earth

if __name__ == '__main__':
    # ISS
    c0 = [cb['radius'] + 414.0, 0.0006189, 51.6393, 0.0, 234.1955, 105.6372]

    # GEO Sattelite
    c1 = [cb['radius'] + 35800, 0, 89, 0, 0, 90]

    # random

    c2 = [cb['radius'] + 35800, .1, 20, 0, 15, 40]

    op0 = OP(c0, tspan, dt, coes=True)
    op1 = OP(c1, tspan, dt, coes=True)
    op2 = OP(c2, tspan, dt, coes=True)

    op0.propogate_orbit()
    op1.propogate_orbit()
    op2.propogate_orbit()

    t.plot_n_orbits([op0.rs, op1.rs, op2.rs], labels=[
                    'ISS', 'GEO', 'random'], show_plot=True, title="Keplerian Orbits")
