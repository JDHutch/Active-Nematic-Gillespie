import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.widgets import Slider
from matplotlib.colors import LinearSegmentedColormap

## Latex Fonts ##
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('legend', fancybox=False, edgecolor='k', shadow=False)
plt.rc('axes', linewidth=3)

fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Report1_0', '1d_u_om_0')

n_surf = np.loadtxt(fold+"/N_Final.dat")

u_unq = np.unique(n_surf[:,0])
q_unq = np.unique(n_surf[:,1])


t_unq = np.unique(n_surf[:,-1])
t_num = t_unq.shape[0]


u_num = u_unq.shape[0]
q_num = q_unq.shape[0]
n_surf[np.where(n_surf[:,2]<0)]
mapper = cm.ScalarMappable(Normalize(np.amin(n_surf[:,2]), np.amax(n_surf[:,2]), clip=False), cm.plasma) # N
#cmap = LinearSegmentedColormap.from_list("", ["royalblue", "crimson"])
#mapper = cm.ScalarMappable(Normalize(np.amin(n_surf[:,2]*(n_surf[:,0]+0.5*n_surf[:,1])), np.amax(n_surf[:,2]*(n_surf[:,0]+0.5*n_surf[:,1])), clip=False), cmap)  # Stress

##### Plot initial surface #####
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax_t = plt.axes([0.95, 0.2, 0.03, 0.65])
t_slider = Slider(ax=ax_t, label=r'$t$', valinit=0, valstep=1, valmin=0, valmax=t_num-1, orientation='vertical')
t_slider.label.set_fontsize(24)

def update(val):
    
    t = t_unq[t_slider.val]

    
    #surf_iter = n_surf[int(t*20**4):int((t+1)*20**4),:]
    surf_iter = n_surf[(n_surf[:,-1] == t),:][1:,:]
    #surf_iter[:,2] *= surf_iter[:,0] + 0.5*surf_iter[:,1]

    sqrLen = int(np.sqrt(surf_iter.shape[0]))
    x = surf_iter[:,0].reshape(sqrLen, sqrLen)
    y = surf_iter[:,1].reshape(sqrLen, sqrLen)
    N = surf_iter[:,2].reshape(sqrLen, sqrLen)

    ax.clear() # Clear axes
    ax.plot_surface(x, y, N, facecolors=mapper.to_rgba(N), lw=0, antialiased=False)
    ax.scatter([1,-1], [1,-1], [np.amin(n_surf[:,2]), np.amax(n_surf[:,2])], alpha=0)
    ax.scatter(x.flatten(), y.flatten(), N.flatten(), c=mapper.to_rgba(N.flatten()))
    ax.set_xlabel(r"$u$", size=18)
    ax.set_ylabel(r"$q$", size=18)
    ax.set_zlabel(r'$n(u,q)$', size=18)


t_slider.on_changed(update)

plt.show()