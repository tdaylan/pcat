import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def _blit_draw(self, artists, bg_cache):
    # Handles blitted drawing, which renders only the artists given instead
    # of the entire figure.
    updated_ax = []
    for a in artists:
        # If we haven't cached the background for this axes object, do
        # so now. This might not always be reliable, but it's an attempt
        # to automate the process.
        if a.axes not in bg_cache:
            # bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            # change here
            bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.figure.bbox)
        a.axes.draw_artist(a)
        updated_ax.append(a.axes)

    # After rendering all the needed artists, blit each axes individually.
    for ax in set(updated_ax):
        # and here
        # ax.figure.canvas.blit(ax.bbox)
        ax.figure.canvas.blit(ax.figure.bbox)

# MONKEY PATCH!!
matplotlib.animation.Animation._blit_draw = _blit_draw

vls = np.linspace(0,2*2*np.pi,100)

fig=plt.figure()
img, = plt.plot(np.sin(vls))
ax = plt.axes()
ax.set_xlim([0,2*2*np.pi])
#ttl = ax.set_title('',animated=True)
ttl = ax.text(.5, 1.05, '', transform = ax.transAxes, va='center')

def init():
    ttl.set_text('')
    img.set_data([0],[0])
    return img, ttl

def func(n):
    ttl.set_text(str(n))
    img.set_data(vls,np.sin(vls+.02*n*2*np.pi))
    return img, ttl

ani = animation.FuncAnimation(fig,func,init_func=init,frames=50,interval=30,blit=True)

plt.show()
