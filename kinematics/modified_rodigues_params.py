import numpy as np
import matplotlib.pyplot as plt

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

def mrp_to_dcm(mrp):

    s1 = mrp[0]
    s2 = mrp[1]
    s3 = mrp[2]

    sn = np.linalg.norm(mrp)

    den = (1 + sn**2)**2
    s_hat = hat(mrp)

    dcm = np.eye(3) + (8 * s_hat @ s_hat - 4 * (1 - sn**2) * s_hat) / den

    return dcm


mrp1 = np.array([0.1,0.2,0.3])

dcm1 = mrp_to_dcm(mrp1)
print(dcm1)

mrp2 = np.array([-0.5,0.1,0.2])

dcm2 = mrp_to_dcm(mrp2)
print(dcm2)

def add_mrp(sv1, sv2):

    sn1 = np.linalg.norm(sv1)
    sn2 = np.linalg.norm(sv2)

    num = (1 - sn2**2) * sv1 + (1 - sn1**2) * sv2 - 2 * np.cross(sv1, sv2)
    den = 1 + sn1**2 * sn2**2 - 2 * np.dot(sn2, sn1)

    return num / den

def sub_mrp(sv1, sv2):

    sn1 = np.linalg.norm(sv1)
    sn2 = np.linalg.norm(sv2)

    num = (1 - sn2**2) * sv1 - (1 - sn1**2) * sv2 + 2 * np.cross(sv1, sv2)
    den = 1 + sn1**2 * sn2**2 + 2 * np.dot(sn2, sn1)

    return num / den

s_B_N = np.array([0.1, 0.2, 0.3])
s_R_B = np.array([-0.1, 0.3, 0.1])
# s_B_R = -s_R_B

s_R_N = add_mrp(s_R_B, s_B_N)
print(s_R_N)

# Problem 3

s_B_N = np.array([0.1, 0.2, 0.3])
s_R_N = np.array([0.5, 0.3, 0.1])

s_B_R = sub_mrp(s_B_N, s_R_N)
print(s_B_R)

#####################################################################################################

def mrp_deriv_mat(sigma, omega):

    s1 = sigma[0]
    s2 = sigma[1]
    s3 = sigma[2]
    sn = np.linalg.norm(sigma)

    sigma_dot = 0.25 * np.array([
        [1 - sn**2 + 2*s1**2,   2*(s1*s2 - s3),         2*(s1*s3 + s2)],
        [2*(s1*s2 + s3),        1 - sn**2 + 2*s2**2,    2*(s2*s3 - s1)],
        [2*(s1*s3 - s2),        2*(s3*s2 + s1),         1 - sn**2 + 2*s3**2]
    ]) @ omega

    return sigma_dot

def mrp_deriv(sigma, omega):

    sn = np.linalg.norm(sigma)
    s_hat = hat(sigma)

    sigma_dot = 0.25 * ((1 - sn**2) * np.eye(3) + 2*s_hat + 2*np.outer(sigma, sigma)) @ omega

    return sigma_dot

# Integrate to get attitude

T = 42.0
dt = 0.001
N = round(T / dt) + 1
t = np.linspace(0.0, T, N)

# omega_t = np.array([np.sin(0.1*t), 0.01*np.ones(N), np.cos(0.1*t)]) * 3.0 * np.pi / 180.0

omega_t = np.zeros((N,3))

for i in range(N):
    omega_t[i,0] = np.sin(0.1*t[i])     * 20.0*np.pi/180.0
    omega_t[i,1] = 0.01                 * 20.0*np.pi/180.0
    omega_t[i,2] = np.cos(0.1*t[i])     * 20.0*np.pi/180.0

# print(omega_t)    

sigma_0 = np.array([0.4, 0.2, -0.1])
sigma_t = np.zeros((N, len(sigma_0)))
sigma_t[0]= sigma_0


for i in range(1, N):

    # calculate derivative of quaternion
    sigma_dot = mrp_deriv_mat(sigma_t[i-1], omega_t[i-1])
    # crp_dot = crp_deriv_mat(crp_t[i-1], omega_t[i-1])

    # integrate by one timestep
    sigma_t[i] = sigma_t[i-1] + sigma_dot * dt

    sigma_norm = np.linalg.norm(sigma_t[i])
    # print(sigma_norm)

    if (np.linalg.norm(sigma_t[i])) > 1.0:
        sigma_t[i] = -sigma_t[i] / sigma_norm

# Print answer to Concept Check 13

sigma_42_norm = np.linalg.norm(sigma_t[N-1])
print(sigma_42_norm)

plt.plot(t, sigma_t)
plt.xlabel("Time (sec)")
plt.ylabel("Modified Rodrigues Parameters")
plt.title("MRP Integration")
plt.legend(["s1", "s2", "s3"])
plt.show()

