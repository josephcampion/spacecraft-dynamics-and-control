# concept_check_9.py

import numpy as np
import matplotlib.pyplot as plt

# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"]})
# for Palatino and other serif fonts use:
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif",
#     "font.serif": ["Palatino"],
# })

def euler_params_deriv(beta, omega):

    b0 = beta[0]
    b1 = beta[1]
    b2 = beta[2]
    b3 = beta[3]

    B_beta = np.array([
        [b0,    -b1,    -b2,    -b3],
        [b1,    b0,     -b3,     b2],
        [b2,    b3,    b0,     -b1],
        [b3,    -b2,     b1,    b0]
    ])

    ep_dot = 0.5 * B_beta @ np.array([0.0,omega[0],omega[1],omega[2]])

    return ep_dot

def euler_params_deriv2(beta, omega):

    w1 = omega[0]
    w2 = omega[1]
    w3 = omega[2]

    W_mat = np.array([
        [0.0,    -w1,    -w2,    -w3],
        [w1,    0.0,     w3,     -w2],
        [w2,    -w3,    0.0,     w1],
        [w3,    w2,     -w1,    0.0]
    ])

    ep_dot = 0.5 * W_mat @ beta

    return ep_dot

beta0 = np.array([0.408248,0.,0.408248,0.816497])

T = 42.0
dt = 0.001
N = round(T / dt) + 1
t = np.linspace(0.0, T, N)

omega_t = np.zeros((N,3))

for i in range(N):
    omega_t[i,0] = np.sin(0.1*t[i])     * 20.0*np.pi/180.0
    omega_t[i,1] = 0.01                 * 20.0*np.pi/180.0
    omega_t[i,2] = np.cos(0.1*t[i])     * 20.0*np.pi/180.0

beta_t = np.zeros((N, len(beta0)))
beta_t[0]= beta0
# print(beta_t)

for i in range(1, N):

    # calculate derivative of quaternion
    ep_dot = euler_params_deriv(beta_t[i-1], omega_t[i-1])

    # integrate by one timestep
    beta_t[i] = beta_t[i-1] + ep_dot * dt


# Print answer to Concept Check 8:
print(beta_t[N-1])

beta_42 = beta_t[N-1]

norm_b1_b2_b3 = np.linalg.norm([beta_42[1], beta_42[2], beta_42[3]])
print(norm_b1_b2_b3)


# plt.plot(t, omega_t)
# plt.xlabel("Time (sec)")
# plt.ylabel("Angular Velocity (rad/s)")
# plt.title("Euler Params Integration")
# plt.legend(["w_x", "w_y", "w_z"])

plt.plot(t, beta_t)
plt.xlabel("Time (sec)")
plt.ylabel("Euler Parameters")
plt.title("Euler Params Integration")
plt.legend(["b0", "b1", "b2", "b3"])
plt.show()

