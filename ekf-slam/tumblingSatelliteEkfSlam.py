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

# class SatelliteSimulation():

#     def __init__(self, satellite, period, dt,)

Jx = 10
Jy = 50
Jz = 50

J = np.array([Jx, Jy, Jz])
# print(J)

T = 10
dt = 1e-2
n = int(T / dt) + 1
# print(n)

t = np.linspace(0.0, T, n)

# print(t)

sat_quat = np.zeros((n,4))
# print(sat_quat)

sat_ang_vel = np.zeros((n,3))

sat_tau = np.array([np.ones(n), np.cos(t), -2*np.sin(t)])

# print(sat_tau)


# plt.plot(t, sat_tau[0])
# plt.plot(t, sat_tau[1])
# plt.plot(t, sat_tau[2])
# plt.show()

sat = Satellite(J, omega_0=[1., 2. ,3.])
print(sat.get_omega_skew())
sat.update_attitude(dt)

# print(sat)
print(sat.get_attitude_quat())

for i in range(n):
    # print(sat_quat[1])
    # sat_quat[i] = sat.get_attitude_quat()
    # print(sat.get_ang_vel())
    sat_ang_vel[i] = sat.get_ang_vel()
    sat_quat[i] = sat.get_attitude_quat()

    sat.update_dynamics(sat_tau[0][i], sat_tau[1][i], sat_tau[2][i], dt)
    sat.update_attitude(dt)

plt.plot(t, sat_quat)
# plt.plot(t, sat_ang_vel)
plt.show()




