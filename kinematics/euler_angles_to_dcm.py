# Euler Angle Addition and Subtraction


import numpy as np


def dir_cosine_mat(q1, q2, q3, n1, n2, n3):
    
    Ri = dcm_helper(q1, n1)
    Rj = dcm_helper(q2, n2)
    Rk = dcm_helper(q3, n3)

    DCM = Rk @ Rj @ Ri

    return DCM

def my_rot_x(q):

    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(q), np.sin(q)],
        [0, -np.sin(q), np.cos(q)]
    ])

    return Rx

def my_rot_y(q):

    Ry = np.array([
        [np.cos(q), 0, -np.sin(q)],
        [0, 1, 0],
        [np.sin(q), 0, np.cos(q)]
    ])

    return Ry

def my_rot_z(q):

    Rz = np.array([
        [np.cos(q), np.sin(q), 0],
        [-np.sin(q), np.cos(q), 0],
        [0, 0, 1]
    ])

    return Rz

def dcm_helper(q, n):

    R = np.eye(3)

    if n == 1:
        
        R = my_rot_x(q)

    elif n == 2:

        R = my_rot_y(q)

    elif n == 3:

        R = my_rot_z(q)

    else:

        print("n should be 1, 2, or 3.")

    return R


# Problem 1

# 3-2-1 rotation matrix (DCM) of Euler angles below:

theta1 =      10 /180*np.pi
theta2 =    20 /180*np.pi
theta3 =      30 /180*np.pi

DCM321 = dir_cosine_mat(theta1, theta2, theta3, 3, 2, 1)

# 3-1-3 rotation matrix (DCM) of Euler angles below:

# psi2 =      38.4812 /180*np.pi;
# theta2 =    -9.84655 /180*np.pi;
# phi2 =      17.4952 /180*np.pi;

psi2 =      40.6423 /180*np.pi;
theta2 =    35.5313 /180*np.pi;
phi2 =      -36.0524 /180*np.pi;
 
DCM_313 = dir_cosine_mat(psi2, theta2, phi2, 3, 1, 3)

# print(DCM_313)

# Problem 2

DCM_BN = dir_cosine_mat(theta1, theta2, theta3, 3, 2, 1)

DCM_RN = dir_cosine_mat(-5*np.pi/180, 5*np.pi/180, 5*np.pi/180, 3, 2, 1) 

DCM_BR = DCM_BN * np.transpose(DCM_RN)

# print(DCM_BN)

# print(DCM_RN)

# print(DCM_BR)

