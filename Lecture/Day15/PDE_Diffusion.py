import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

def double_deriv(dx, C):
    D = 0.5
    dCd2 = np.zeros(len(C))
    for i in range(1, len(C) - 1):
        dCd2[i] = D * (C[i + 1] - 2 * C[i] + C[i - 1]) / (dx * dx)
    return dCd2

def deriv(dx, C):
    v = 0.25
    dC = np.zeros(len(C))
    dC[0] = -1 + 0.5 * C[0]
    dC[1:-1] = v * (C[2:] - C[:-2]) / (2 * dx)
    dC[-1] = 1 + 0.5 * C[-1]
    return dC

def main(dx, L, dt, t_f):
    N_steps = int(t_f / dt)
    N = int(L / dx)
    C = np.ones((N, N_steps))
    for i in range(1, N_steps):
        dC = deriv(dx, C[:, i - 1])
        dC2 = double_deriv(dx, C[:, i - 1])
        C[:, i] = C[:, i - 1] + dt * (0.5 * dC - 0.25 * dC2 - C[:, i - 1])
    return C

L = 100
dx = 1
t_f = 10
dt = 0.1

C = main(dx, L, dt, t_f)
x = np.linspace(0, L, int(L / dx))

metadata = dict(title='Diffusion', artist='Matplotlib')
writer = FFMpegWriter(fps=15, codec='h264', bitrate=3000, metadata=metadata)

fig, ax = plt.subplots()
ax.set_xlim(x[0], x[-1])
ax.set_ylim(0, 1)
line, = ax.plot([], [], 'b')

with writer.saving(fig, "diffusion_simulation.mp4", dpi=100):
    for i in range(0, C.shape[1], 5):
        line.set_data(x, C[:, i])
        ax.set_title(f'Time: {i * dt:.2f} s')
        writer.grab_frame()
