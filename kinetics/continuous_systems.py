import numpy as np

R1 = np.array([1, -1, 2])
R2 = np.array([-1, -3, 2])
R3 = np.array([2, -1, -1])
R4 = np.array([3, -1, -2])

m1 = 1
m2 = 1
m3 = 2
m4 = 2

M = (m1 + m2 + m3 + m4)

R1dot = np.array([2, 1, 1])
R2dot = np.array([0, -1, 1])
R3dot = np.array([3, 2, -1])
R4dot = np.array([0, 0, 1])

Rcm = (m1*R1 + m2*R2 + m3*R3 + m4*R4) / M
print(Rcm)

Rcm_dot = (m1*R1dot + m2*R2dot + m3*R3dot + m4*R4dot) / M
print(Rcm_dot)

r1 = R1 - Rcm
r2 = R2 - Rcm
r3 = R3 - Rcm
r4 = R4 - Rcm

r1dot = R1dot - Rcm_dot
r2dot = R2dot - Rcm_dot
r3dot = R3dot - Rcm_dot
r4dot = R4dot - Rcm_dot

T_trans = 0.5 * M * np.dot(Rcm_dot, Rcm_dot)
print(T_trans)


#################

Hc = m1 * np.cross(r1, r1dot) + \
        m2 * np.cross(r2, r2dot) + \
        m3 * np.cross(r3, r3dot) + \
        m4 * np.cross(r4, r4dot)
print(Hc)

H0 = M * np.cross(Rcm, Rcm_dot) + Hc
print(H0)

