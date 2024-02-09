import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
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

if __name__ == "__main__":
    perts = null_perts()
    perts['J2'] = True
    op = OP(t.tle2coes('iss.txt'), tspan, dt, perts=perts, coes=True)

    # t.plot_n_orbits([op.rs], show_plot=True, labels=['perts'])

    op.calculate_coes()