import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font_size = 15
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 500,
    'axes.labelsize': font_size,
    'axes.titlesize': font_size,
    'lines.markersize': 3,
    'lines.linewidth': 1.4,
    'font.size': font_size,
    'legend.fontsize': 13,
    'xtick.labelsize': font_size,
    'text.usetex': True,
    'ytick.labelsize': font_size,
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)

sample_path = os.path.join('singleGraph', '3', 'line_T.xy')
cases = [
    dict(
        id='explicit',
        style=dict(label= 'expl. RK4', color='k' , linestyle='solid'),
    ),
    dict(
        id='implicit-Euler',
        style=dict(label='impl. Euler', color='tab:red', linestyle='dashed'),
    ),
    dict(
        id='implicit-backward',
        style=dict(label='impl. backward', color='tab:purple'  , linestyle='dashdot'),
    ),
    dict(
        id='implicit-Crank-Nicolson',
        style=dict(label='impl. C-N 0.9', color='tab:blue'  , linestyle='dotted'),
    ),

]

plt.figure()

# analytical solution
L = 2*np.pi
x = np.linspace(-L/2, L/2, 200)
sigma = L/20
x0 = 0
T = lambda x: np.exp(-((x-x0)*(x-x0))/(2*sigma*sigma))
plt.plot(x, T(x), color='tab:gray', linestyle='dashed', label='analytical')

for case in cases:
    data_path = os.path.join('postProcessing-' + case['id'], sample_path)
    d = np.loadtxt(data_path, unpack=True)
    plt.plot(d[0], d[1], **case['style'])


plt.xlabel('x, m')
plt.ylabel('T, [-]')
plt.xticks(ticks=[-np.pi, -0.5*np.pi, 0, 0.5*np.pi, np.pi], labels=['$-\pi$', '$-\pi/2$', '0', '$\pi$/2', '$\pi$'])
plt.legend()
plt.grid()
plt.savefig('convection1d_ddt.pdf')
