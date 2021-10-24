import sys
sys.path.append('..')
import numpy as np

from kinematics.euler_angles_to_dcm import dir_cosine_mat
from kinematics.modified_rodrigues_params import hat

I_C_B = np.array([
    [10, 1, -1],
    [1, 5, 1],
    [-1, 1, 8]
])

m = 12.5

phi = -10 * np.pi / 180
theta = 10 * np.pi / 180
psi = 5 * np.pi / 180

RcP_N = np.array([-0.5, 0.5, 0.25])

R_B_N = dir_cosine_mat(phi, theta, psi, 3, 2, 1)

Rc_hat = hat(R_B_N @ RcP_N)

I_P_B = I_C_B + m * Rc_hat @ Rc_hat.transpose()

print(I_P_B)


