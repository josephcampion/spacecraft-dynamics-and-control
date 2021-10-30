import numpy as np 
import matplotlib.pyplot as plt 

# What is the definiteness of the following matrices?

K1 = np.array([
    [1.53947, -0.0422688, -0.190629],
    [-0.0422688, 1.4759, 0.459006],
    [-0.190629, 0.459006, 1.48463]
])

w1, v1 = np.linalg.eig(K1)

print(w1)
print(v1)
print()

K2 = np.array([
    [-0.984331, -1.10006, -0.478579],
    [-1.10006, 1.03255, 0.338318],
    [-0.478579, 0.338318, 1.45178]
])

w2, v2 = np.linalg.eig(K2)

print(w2)
print(v2)

K3 = np.array([
    [-2.0353, 0.296916, -0.365128],
    [0.296916, -1.10369, -0.074481],
    [-0.365128, -0.074481, -2.86101]
])

w3, v3 = np.linalg.eig(K3)

print(w3)
print(v3)

# Plot V(x) = ln(1 + x1**2 + x2**2)

n = 100
x = np.linspace(-10, 10, n)
y = np.linspace(-10, 10, n)
xv, yv = np.meshgrid(x, y)

Vx = np.log(np.ones((n,n)) + xv**2 + yv**2)

fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot3D(xv, yv, Vx)
plt.contourf(xv, yv, Vx)
plt.show()

