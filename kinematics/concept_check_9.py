# concept_check_9.py

import numpy as np
import matplotlib.pyplot as plt

def euler_deriv(theta, phi, omega):

    mat_A = np.array([
            [0.0, np.sin(phi), np.cos(phi)],
            [0.0, np.cos(phi) * np.cos(theta), -np.sin(phi) * np.cos(theta)],
            [np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi) * np.sin(theta)]
        ])

    eul_dot = (1 / np.cos(theta)) * mat_A @ omega

    return eul_dot

psi0 = 40.0 * np.pi / 180
theta0 = 30.0 * np.pi / 180
phi0 = 80.0 * np.pi / 180

x0 = np.array([psi0, theta0, phi0])

T = 42.0
dt = 0.1
N = round(T / dt) + 1
t = np.linspace(0.0, T, N)

omega_t = np.array([np.sin(0.1*t), 0.01*np.ones(N), np.cos(0.1*t)]) * 20 * np.pi / 180

xt = np.zeros((N, len(x0)))

xt[0]= x0

for i in range(2, N):

    theta_i_minus_1 = xt[i-1,1]
    phi_i_minus_1 = xt[i-1,2]

    eul_dot = euler_deriv(theta_i_minus_1, phi_i_minus_1, omega_t[0:3,i-1])
    
    xt[i] = xt[i-1] + np.transpose(eul_dot) * dt

# plt.plot(t, np.transpose(omega_t))
plt.plot(t,xt)
plt.show()





