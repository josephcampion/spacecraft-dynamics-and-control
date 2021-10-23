import numpy as np
from concept_check_5_6 import dcm_to_eul_params, eul_params_to_dcm

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