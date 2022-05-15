import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


def set_u0(x):
    if 0.4 <= x < 0.6:
        return 4
    else:
        return 1


N = 200
T = 0.01

a = 1
b = 0
c = 1

h = (a - b) / N

dt = h*h / (2 * c) * 0.9

space = np.linspace(a, b, N)

u0 = np.array(list(map(set_u0, space)))
u1 = np.zeros(u0.shape)

# print(u0)
# plt.plot(u0)

boundary = 'Neumann'

t = 0
while t < T:
    for i in range(1, N - 1):
        u1[i] = u0[i] + c * (dt / (h*h)) * (u0[i-1] - 2*u0[i] + u0[i+1])
        if boundary == 'Dirichlet':
            u1[0] = u0[0]
            u1[N-1] = u0[N-1]
        elif boundary == 'Neumann':
            u1[0] = u0[1]
            u1[N - 1] = u0[N - 2]

    u0 = u1
    t += dt


plt.title(f" Homogeneous {boundary}. N={N}, T={T}")
plt.plot(u0)
plt.show()
