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

def quat_to_skew_mat(beta):

    b0 = beta[0]
    b1 = beta[1]
    b2 = beta[2]
    b3 = beta[3]

    B = np.array([
        [b0,    -b1,    -b2,    -b3],
        [b1,    b0,     b3,     -b2],
        [b2,    -b3,    b0,     b1],
        [b3,    b2,     -b1,    b0]
    ])

    return B


beta_BN = np.array([0.774597,0.258199,0.516398,0.258199])
beta_FB = np.array([0.359211,0.898027,0.179605,0.179605])

mat_FB = quat_to_skew_mat(beta_FB)
print(mat_FB)

beta_FN = mat_FB @ beta_BN

print(beta_FN)


####################################33

beta_FN = np.array([0.359211,0.898027,0.179605,0.179605])
beta_BN = np.array([-0.377964,0.755929,0.377964,0.377964])

dcm_FN = eul_params_to_dcm(beta_FN)
dcm_BN = eul_params_to_dcm(beta_BN)

dcm_NB = np.transpose(dcm_BN)

dcm_FB = dcm_FN @ dcm_NB

beta_FB = dcm_to_eul_params(dcm_FB)

print(beta_FB)

#######################################

# TODO: fix below...
# beta_NB = np.array([-0.377964,-0.755929,-0.377964,-0.377964])

# mat_FN = quat_to_skew_mat(beta_FN)

# beta_FB = mat_FN @ beta_NB

# print(beta_FB)