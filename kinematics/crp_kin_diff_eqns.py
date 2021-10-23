import numpy as np
import matplotlib.pyplot as plt

# CRP differential kinematic equations

def hat(omega):
    """
    vector to skew symmetric matrix associated with cross product
    """
    a = omega.item(0)
    b = omega.item(1)
    c = omega.item(2)

    omega_hat = np.array([[0, -c, b],
                          [c, 0, -a],
                          [-b, a, 0]])
    return omega_hat

def crp_deriv(crp, omega):

    q1 = crp[0]
    q2 = crp[1]
    q3 = crp[2]

    w1 = omega[0]
    w2 = omega[1]
    w3 = omega[2]

    crp_dot = 0.5 * np.array([
        [1 + q1**2,     q1*q2 - q3,     q1*q3 + q2],
        [q2*q1 + q3,    1 + q2**2,      q2*q3 - q1],
        [q3*q1 - q2,    q3*q2 + q1,     1 + q3**2]
    ]) @ omega

    return crp_dot

def crp_deriv_mat(crp, omega):

    q_hat = hat(crp)

    crp_dot = 0.5 * (np.eye(3) + q_hat + np.outer(crp,crp)) @ omega

    return crp_dot

# Integrate to get attitude

T = 42.0
dt = 0.01
N = round(T / dt) + 1
t = np.linspace(0.0, T, N)

# omega_t = np.array([np.sin(0.1*t), 0.01*np.ones(N), np.cos(0.1*t)]) * 3.0 * np.pi / 180.0

omega_t = np.zeros((N,3))

for i in range(N):
    omega_t[i,0] = np.sin(0.1*t[i])     * 3.0*np.pi/180.0
    omega_t[i,1] = 0.01                 * 3.0*np.pi/180.0
    omega_t[i,2] = np.cos(0.1*t[i])     * 3.0*np.pi/180.0

# print(omega_t)    

crp0 = np.array([0.4, 0.2, -0.1])
crp_t = np.zeros((N, len(crp0)))
crp_t[0]= crp0


# print(crp_t)


for i in range(1, N):

    # calculate derivative of quaternion
    crp_dot = crp_deriv(crp_t[i-1], omega_t[i-1])
    # crp_dot = crp_deriv_mat(crp_t[i-1], omega_t[i-1])

    # integrate by one timestep
    crp_t[i] = crp_t[i-1] + crp_dot * dt

# Print answer to Concept Check 13

q_42_norm = np.linalg.norm(crp_t[N-1])
print(q_42_norm)


# Plot results

# plt.plot(t, omega_t)
# plt.xlabel("Time (sec)")
# plt.ylabel("Angular Velocity (rad/s)")
# plt.title("Euler Params Integration")
# plt.legend(["w_x", "w_y", "w_z"])
# plt.show()

plt.plot(t, crp_t)
plt.xlabel("Time (sec)")
plt.ylabel("Classic Rodrigues Parameters")
plt.title("CRP Integration")
plt.legend(["q0", "q1", "q2"])
plt.show()



