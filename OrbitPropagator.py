import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import OrbitTools as t

import planetary_data as pd


def null_perts():
    return {
        'J2': False,
        'aero': False,
        'moon_grav': False,
        'solar_gravity': False
    }


class OrbitPropagator:
    def __init__(self, state0, tspan, dt, coes=False, deg=True, cb=pd.earth, perts=null_perts()):
        if coes:
            self.r0, self.v0 = t.coes2rv(state0, deg=True, mu=cb['mu'])
        else:
            self.r0 = state0[:3]

            self.v0 = state0[3:]

        self.y0 = self.r0.tolist() + self.v0.tolist()
        self.tspan = tspan
        self.dt = dt
        self.cb = cb

        self.n_steps = int(np.ceil(self.tspan/self.dt))

        self.ys = np.zeros((self.n_steps, 6))
        self.ts = np.zeros((self.n_steps, 1))

        self.ys[0] = np.array(self.y0)

        self.step = 1

        self.solver = ode(self.diffy_q)
        self.solver.set_integrator('lsoda')

        self.solver.set_initial_value(self.y0, 0)

        self.perts = perts

        self.propogate_orbit()

        # self.solver.set_f_params(self.cb['mu']) idk why you dont have to set the mu value

    def propogate_orbit(self):

        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t + self.dt)
            self.ts[self.step] = self.solver.t
            self.ys[self.step] = self.solver.y
            self.step += 1

        self.rs = self.ys[:, :3]
        self.vs = self.ys[:, 3:]

    def diffy_q(self, t, y):
        # unpack state

        rx, ry, rz, vx, vy, vz = y

        r = np.array([rx, ry, rz])
        v = np.array([vx, vy, vz])
        norm_r = np.linalg.norm(r)

        a = -r*self.cb['mu']/norm_r**3

        # J2 Pertubation
        if self.perts['J2']:
            z2 = r[2]**2
            r2 = norm_r**2
            tx = r[0]/norm_r*(5*z2/r2-1)
            ty = r[1]/norm_r*(5*z2/r2-1)
            tz = r[2]/norm_r*(5*z2/r2-3)

            a_j2 = 1.5*self.cb['J2']*self.cb['mu'] * \
                self.cb['radius']**2/norm_r**4*np.array([tx, ty, tz])

        return [vx, vy, vz, a[0], a[1], a[2]]

    def plot_3d(self, show_plot=False, save_plot=False, title='No Title Given'):
        fig = plt.figure(figsize=(16, 8))
        ax = fig.add_subplot(111, projection='3d')

        ax.plot(self.rs[:, 0], self.rs[:, 1], self.rs[:, 2], 'b', label="traj")
        ax.plot([self.rs[0, 0]], [self.rs[0, 1]], [
                self.rs[0, 2]], 'ro', label="initial pos")

        _u, _v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        _x = self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y = self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z = self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x, _y, _z, cmap="Blues")

        l = self.cb['radius'] * 2

        x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        u, v, w = [[l, 0, 0], [0, l, 0], [0, 0, l]]

        ax.quiver(x, y, x, u, v, w, color="k")

        max_val = np.max(np.abs(self.rs))

        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])

        ax.set_xlabel(['X (km)'])
        ax.set_ylabel(['Y (km)'])
        ax.set_zlabel(['Z (km)'])

        ax.set_aspect('auto')

        ax.set_title(title)

        plt.legend()
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title + '.png', dpi=300)
