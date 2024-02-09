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
    op0 = OP(t.tle2coes('iss.txt'), tspan, dt, coes=True, deg=False)
    op1 = OP(t.tle2coes('hst.txt'), tspan, dt, coes=True, deg=False)

    t.plot_n_orbits([op0.rs, op1.rs], labels=['iss', 'hst'], show_plot=True)
