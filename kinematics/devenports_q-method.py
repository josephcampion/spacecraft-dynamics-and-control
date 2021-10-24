import numpy as np
from concept_check_5_6 import eul_params_to_dcm

def dp_q_method_2d(vm1_B, vm2_B, vm1_N, vm2_N):

    w1 = 1.0
    w2 = 1.0

    B = w1 * np.outer(vm1_B, vm1_N) + w2 * np.outer(vm2_B, vm2_N)

    S = B + B.transpose()

    sigma = np.trace(B)

    Z = np.array([
        B[1,2] - B[2,1],
        B[2,0] - B[0,2],
        B[0,1] - B[1,0]
    ])

    K = np.zeros((4,4))
    K[0,:]= np.concatenate((sigma, Z), axis=None)
    K[1:4,0] = Z
    K[1:4,1:4] = S - sigma * np.eye(3)

    w,v = np.linalg.eig(K)
    print(w)
    print(v)

    eig_max_ind = np.argmax(w)
    beta = v[:,eig_max_ind]
    print(beta)

    return eul_params_to_dcm(beta)

# Problem 1

vm1_B = np.array([0.8273,       0.5541,     -0.0920])
vm2_B = np.array([-0.8285,      0.5522,     -0.0955])
vm1_N = np.array([-0.1517,      -0.9669,    0.2050])
vm2_N = np.array([-0.8393,      0.4494,     -0.3044])

BN = dp_q_method_2d(vm1_B, vm2_B, vm1_N, vm2_N)
print(BN)

