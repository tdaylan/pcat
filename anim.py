import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

figr, axis = plt.subplots()

axis.set_xlim([0., 10.])
axis.set_ylim([0., 10.])
axis.set_aspect('equal')

patch = plt.Circle((5, -5), 1., color='y')

def init():
    patch.center = (5, 5)
    axis.add_patch(patch)
    return patch,

def animate(i):
    x = 5 + 3 * np.sin(np.radians(i))
    y = 5 + 3 * np.cos(np.radians(i))
    patch.center = (x, y)
    return patch,

anim = animation.FuncAnimation(figr, animate, 
                               init_func=init, 
                               frames=360, 
                               interval=20,
                               blit=True)

plt.show()
