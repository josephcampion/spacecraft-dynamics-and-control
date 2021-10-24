import numpy as np 

def triad_method(vm1_B, vm2_B, vm1_N, vm2_N):

    tv1_B = vm1_B
    vm1_x_vm2_B = np.cross(vm1_B, vm2_B)
    tv2_B =  vm1_x_vm2_B / np.linalg.norm(vm1_x_vm2_B)
    tv3_B = np.cross(tv1_B, tv2_B)

    tv1_N = vm1_N
    vm1_x_vm2_N = np.cross(vm1_N, vm2_N)
    tv2_N = vm1_x_vm2_N / np.linalg.norm(vm1_x_vm2_N)
    tv3_N = np.cross(tv1_N, tv2_N)

    BT = np.array([tv1_B, tv2_B, tv3_B]).transpose()
    # print(BT)

    NT = np.array([tv1_N, tv2_N, tv3_N]).transpose()

    BN = BT @ np.transpose(NT)

    return BN

# Problem 1

vm1_B = np.array([0.8273,       0.5541,     -0.0920])
vm2_B = np.array([-0.8285,      0.5522,     -0.0955])
vm1_N = np.array([-0.1517,      -0.9669,    0.2050])
vm2_N = np.array([-0.8393,      0.4494,     -0.3044])

BN = triad_method(vm1_B, vm2_B, vm1_N, vm2_N)
print(BN)

# Problem 2

BN_est = np.array([
    [0.969846,-0.200706,-0.138258],
    [0.17101,0.96461,-0.200706],
    [0.173648,0.17101,0.969846]
])

BN_true = np.array([
    [0.963592,-0.223042,-0.147454],
    [0.187303,0.956645,-0.223042],
    [0.190809,0.187303,0.963592],
])

BB_est = BN_est @ np.transpose(BN_true)
print(BB_est)

cos_phi = 0.5 * (np.trace(BB_est) - 1)

phi_deg = np.arccos(cos_phi) * 180.0 / np.pi
print(phi_deg)
