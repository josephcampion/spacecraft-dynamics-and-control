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



