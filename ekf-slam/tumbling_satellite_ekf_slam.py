# tumblingSatelliteEkfSlam.py

import numpy as np
import matplotlib.pyplot as plt
# import scipy
# import random
from tools.rotations import Quaternion2Rotation 


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

    def update_attitude(self, dt):
        At = np.array([
            [0., -self.omega[0], -self.omega[1], -self.omega[2]],
            [self.omega[0], 0., -self.omega[2], self.omega[1]],
            [self.omega[1], self.omega[2], 0., -self.omega[0]],
            [self.omega[2], -self.omega[1], self.omega[0], 0.]
        ])
        new_quat = (np.identity(4) + At * dt) @ self.quat
        self.quat = new_quat / np.linalg.norm(new_quat)        

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
        self.sat_omega = np.zeros((n,3))

    def set_attitude(self, attitude_quat):
        self.sat_quat = attitude_quat

    def get_sat_attitude(self):
        return self.sat.sat_quat

    def set_sat_ang_vel(self, new_sat_ang_vel):
        self.sat_omega = new_sat_ang_vel

    def get_sat_ang_vel(self):
        return self.sat_omega

    def set_control(self, new_control):
        self.tau = new_control

    def run_simulation(self, plot_results=True):
        
        # Iterate through each timestep
        for i in range(self.n):

            # store attitude and angular velocity
            self.sat_omega[i] = self.sat.get_ang_vel()
            self.sat_quat[i] = self.sat.get_attitude_quat()

            # update dynamics and kinematics based on control effort
            self.sat.update_dynamics(self.tau[0][i], self.tau[1][i], self.tau[2][i], dt)
            self.sat.update_attitude(dt)
        
        if (plot_results):
            fig, (ax1, ax2) = plt.subplots(2,1)
            for i in range(4):
                ax1.plot(self.t, self.sat_quat[0:, i], label=("e"+str(i)))
            for i in range(3):
                ax2.plot(self.t, self.sat_omega[0:, i], label=("w"+str(i+1)))
            fig.suptitle('Satellite Attitude')
            ax1.set_ylabel('Quaternion Value')
            ax1.legend()
            ax2.set_ylabel('Omega [rad/s]')
            ax2.set_xlabel('Time [s]')
            ax2.legend()
            plt.show()

################################################################################################

class StarSensor():

    def __init__(self, number_of_stars=10, star_coords=None, meas_noise=1e-4):
        
        if np.any(star_coords):
            azimuth_angles = star_coords[:,0]
            polar_angles = star_coords[:,1]
            self.N = np.size(star_coords,0)
            if number_of_stars != self.N:
                print("Using size of star coordinates to determine number of stars.")
        else:
            self.N = number_of_stars
            azimuth_angles = 2 * np.pi * np.random.random_sample((self.N,1))
            polar_angles = np.pi * np.random.random_sample((self.N,1))
            
        self.phi = azimuth_angles
        self.theta = polar_angles
        self.vt = meas_noise

    def set_azimuthal_angles(self, new_azimuthal_angles):
        self.phi = new_azimuthal_angles

    def get_azimuthal_angles(self):
        return self.phi 

    def set_polar_angles(self, new_polar_angles):
        self.theta = new_polar_angles

    def get_polar_angles(self):
        self.theta

    def get_unit_vectors(self):
        sv = np.zeros((self.N,3))
        for i in range(self.N):
            sv[i,0] = np.sin(self.theta[i]) * np.cos(self.phi[i])
            sv[i,1] = np.sin(self.theta[i]) * np.sin(self.phi[i])
            sv[i,2] = np.cos(self.theta[i])
        return sv

################################################################################################

class SatExtendedKalmanFilter():

    def __init__(self, satellite, star_sensor, process_noise, meas_noise):
        self.sat = satellite
        self.ss = star_sensor
        self.update_model_jacob()
        # self.init_jacobians()
        self.Qt = process_noise
        self.Rt = meas_noise

    def update_model_jacob(self):
        omega = self.sat.get_ang_vel()
        Wt = np.array([
            [0., -omega[0], -omega[1], -omega[2]],
            [omega[0], 0., -omega[2], omega[1]],
            [omega[1], omega[2], 0., -omega[0]],
            [omega[2], -omega[1], omega[0], 0.]
        ])
        At = (np.identity(4) + Wt * dt)
        self.At = At

    def update_meas_jacob(self):
        pass

    def get_dg_dq0(self, star_vector):
        sx = star_vector[0]
        sy = star_vector[1]
        sz = star_vector[2]
        q = self.sat.get_attitude_quat()
        # q0 = q[0]
        q1 = q[1]
        q2 = q[2]
        q3 = q[3]
        dg_dq0 = np.array([
            [ 2 * q3 * sy - 2 * q2 * sz ],
            [-2 * q3 * sx + 2 * q1 * sz ],
            [ 2 * q2 * sx - 2 * q1 * sy ]
        ])
        return dg_dq0

    def get_dg_dq1(self, star_vector):
        sx = star_vector[0]
        sy = star_vector[1]
        sz = star_vector[2]
        q = self.sat.get_attitude_quat()
        q0 = q[0]
        q1 = q[1]
        q2 = q[2]
        q3 = q[3]
        dg_dq1 = np.array([
            [ 2 * q2 * sy - 2 * q3 * sz ],
            [ 2 * q2 * sx - 4 * q1 * sy + 2 * q0 * sz ],
            [ 2 * q3 * sx - 2 * q0 * sy - 4 * q1 * sz ]
        ])
        return dg_dq1

    def get_dg_dq2(self, star_vector):
        sx = star_vector[0]
        sy = star_vector[1]
        sz = star_vector[2]
        q = self.sat.get_attitude_quat()
        q0 = q[0]
        q1 = q[1]
        q2 = q[2]
        q3 = q[3]
        dg_dq2 = np.array([
            [-4 * q2 * sx + 2 * q1 * sy + 2 * q0 * sz ],
            [ 2 * q1 * sx + 2 * q3 * sz ],
            [ 2 * q0 * sx + 2 * q3 * sy - 4 * q2 * sz ]
        ])
        return dg_dq2

    def get_dg_dq3(self, star_vector):
        sx = star_vector[0]
        sy = star_vector[1]
        sz = star_vector[2]
        q = self.sat.get_attitude_quat()
        q0 = q[0]
        q1 = q[1]
        q2 = q[2]
        q3 = q[3]
        dg_dq3 = np.array([
            [-4 * q3 * sx + 2 * q0 * sy + 2 * q1 * sz ],
            [-2 * q0 * sx - 4 * q3 * sy + 2 * q2 * sz ],
            [ 2 * q1 * sx + 2 * q2 * sy ]
        ])
        return dg_dq3

################################################################################################

# satellite and simulation params
Jx = 10
Jy = 50
Jz = 50
J = np.array([Jx, Jy, Jz])
T = 10
dt = 1e-2
n = int(T / dt) + 1
t = np.linspace(0.0, T, n)
sat_tau = np.array([np.ones(n), np.cos(t), -2*np.sin(t)])

sat = Satellite(J, omega_0=[1., 2. , 3.])
sat_sim = SatelliteSimulation(sat, T, dt, sat_tau)
sat_sim.run_simulation(False)

ss = StarSensor()

sat_ekf = SatExtendedKalmanFilter(sat, ss, 1e-4, 1e-4)
# print(sat_ekf.At)

sv = np.array([0., np.sqrt(3)/2, 0.5])
print(sat_ekf.get_dg_dq3(sv))