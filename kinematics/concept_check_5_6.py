import numpy as np


def eul_params_to_dcm(beta):

    b0 = beta[0]
    b1 = beta[1]
    b2 = beta[2]
    b3 = beta[3]

    dcm = np.array([
        [b0**2 + b1**2 - b2**2 - b3**2, 2*(b1*b2 + b0*b3),              2*(b1*b3 - b0*b2)],
        [2*(b1*b2 - b0*b3),             b0**2 + b1**2 - b2**2 - b3**2,  2*(b2*b3 + b0*b1)],
        [2*(b1*b3 + b0*b2),             2*(b2*b3 - b0*b1),              b0**2 - b1**2 - b2**2 + b3**2]
    ])

    return dcm

beta = np.array([0.235702,0.471405,-0.471405,0.707107])

dcm = eul_params_to_dcm(beta)

print(dcm)

beta_inv = np.array([0.235702,-0.471405,0.471405,-0.707107])

dcm_inv = eul_params_to_dcm(beta_inv)

print(dcm_inv)

#################################################################3

def dcm_to_eul_params(dcm):

    b0 = 0.5 * np.sqrt(dcm[0,0] + dcm[1,1] + dcm[2,2] + 1)
    b1 = (dcm[1,2] - dcm[2,1]) / (4 * b0)
    b2 = (dcm[2,0] - dcm[0,2]) / (4 * b0)
    b3 = (dcm[0,1] - dcm[1,0]) / (4 * b0)

    beta = np.array([b0, b1, b2, b3])

    return beta


dcm = np.array([
    [-0.529403,-0.474115,0.703525], 
    [-0.467056, -0.529403,-0.708231],
    [0.708231, -0.703525, 0.0588291]
])

beta = dcm_to_eul_params(dcm)

print(beta)

##################################################################

def eul_to_quat(phi, theta, psi):
    """
    Converts an euler angle attitude to a quaternian attitude
    :param euler: Euler angle attitude in a np.matrix(phi, theta, psi)
    :return: Quaternian attitude in np.array(e0, e1, e2, e3)
    """

    e0 = np.cos(psi/2.0) * np.cos(theta/2.0) * np.cos(phi/2.0) + np.sin(psi/2.0) * np.sin(theta/2.0) * np.sin(phi/2.0)
    e1 = np.cos(psi/2.0) * np.cos(theta/2.0) * np.sin(phi/2.0) - np.sin(psi/2.0) * np.sin(theta/2.0) * np.cos(phi/2.0)
    e2 = np.cos(psi/2.0) * np.sin(theta/2.0) * np.cos(phi/2.0) + np.sin(psi/2.0) * np.cos(theta/2.0) * np.sin(phi/2.0)
    e3 = np.sin(psi/2.0) * np.cos(theta/2.0) * np.cos(phi/2.0) - np.cos(psi/2.0) * np.sin(theta/2.0) * np.sin(phi/2.0)

    return np.array([e0, e1, e2, e3])


phi = 20.0 * np.pi / 180.0
theta = 10.0 * np.pi / 180.0
psi = -10.0 * np.pi / 180.0

print(eul_to_quat(psi, theta, phi))
