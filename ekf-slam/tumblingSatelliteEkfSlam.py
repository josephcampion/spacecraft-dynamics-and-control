# tumblingSatelliteEkfSlam.py

import numpy as np
import matplotlib.pyplot as plt
# import scipy


class Satellite():

    def __init__(self, J, quat_0=np.array([1.0, 0.0, 0.0, 0.0]), omega_0=np.array([0.0, 0.0, 0.0])):
        self.Jx = J[0]
        self.Jy = J[1]
        self.Jz = J[2]
        self.quat = quat_0
        self.omega = omega_0

    def get_attitude_quat(self):
        return self.quat

    def get_ang_vel(self):
        return self.omega

    def update_dynamics(self, tau_x, tau_y, tau_z, dt):
        self.omega[0] = self.omega[0] + ((self.Jy - self.Jz) * self.omega[1] * self.omega[2] + tau_x) / self.Jx * dt # + noise
        self.omega[1] = self.omega[1] + ((self.Jz - self.Jx) * self.omega[2] * self.omega[0] + tau_y) / self.Jy * dt # + noise
        self.omega[2] = self.omega[2] + ((self.Jx - self.Jy) * self.omega[0] * self.omega[1] + tau_z) / self.Jz * dt # + noise

    def get_omega_skew(self):
        Omega = np.array([
            [0., -self.omega[0], self.omega[1]],
            [self.omega[0], 0., -self.omega[2]],
            [-self.omega[1], self.omega[2], 0.]
        ])
        return Omega

    def update_attitude(self, dt):
        At = np.array([
            [0., -self.omega[0], -self.omega[1], -self.omega[2]],
            [self.omega[0], 0., -self.omega[0], self.omega[1]],
            [self.omega[1], self.omega[0], 0., -self.omega[2]],
            [self.omega[2], -self.omega[1], self.omega[2], 0.]
        ])
        # print(At)
        Bt = (np.identity(4) + At * dt)
        # print(Bt)
        new_quat = (np.identity(4) + At * dt) @ self.quat
        self.quat = new_quat / np.linalg.norm(new_quat)
        # print(self.quat)
        

    def __str__(self):
        return "Satellite with inertia: [%f, %f, %f]" % ( self.Jx, self.Jy, self.Jz )

##################################################################################################3

class SatelliteSimulation():

    def __init__(self, satellite, period, timestep, control):
        self.sat = satellite
        self.T = period
        self.dt = timestep
        # vector of control effort (torque about x,y,z):
        self.tau = control
        # number of timesteps:
        self.n = int(period / dt) + 1
        self.t = np.linspace(0.0, period, self.n)
        # satellite attitude (quaternion) at each timestep:
        self.sat_quat = np.zeros((n,4))
        # satellite angular velocity (omega) at each timestep:
        self.sat_ang_vel = np.zeros((n,3))

    def run_simulation(self, plot_results=True):
        
        # Iterate through each timestep
        for i in range(self.n):

            # store attitude and angular velocity
            self.sat_ang_vel[i] = self.sat.get_ang_vel()
            self.sat_quat[i] = self.sat.get_attitude_quat()

            # update dynamics and kinematics based on control effort
            self.sat.update_dynamics(self.tau[0][i], self.tau[1][i], self.tau[2][i], dt)
            self.sat.update_attitude(dt)
        
        if (plot_results):
            fig, ax = plt.subplots() #(2,1,1)
            ax.plot(self.t, self.sat_quat)
            # ax.title('Satellite Attitude')
            # ax.ylabel('Quaternion Value')
            plt.show()

            # legend = ax.legend()


################################################################################################

# satellite and simulation params
Jx = 10
Jy = 50
Jz = 50
J = np.array([Jx, Jy, Jz])
T = 10
dt = 1e-2
sat_tau = np.array([np.ones(n), np.cos(t), -2*np.sin(t)])

sat = Satellite(J, omega_0=[1., 2. ,3.])
sat_sim = SatelliteSimulation(sat, T, dt, sat_tau)
sat_sim.run_simulation()


