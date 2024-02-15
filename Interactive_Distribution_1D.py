import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.widgets import Slider, RadioButtons

## Latex Fonts ##
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('legend', fancybox=False, edgecolor='k', shadow=False)
plt.rc('axes', linewidth=3)

def dbb(u, q, mu, D, chi, e0, rho):
    xx = mu*u**2+2*D*u*q+chi*q**2
    return np.exp(-e0-0.5*(xx)/(rho*1E18*300*1.381E-23))

def omega(u, q, a, b, c, mu, D, chi):
    return a*np.exp(-0.5*b*(mu*((u-c)**2)+2*D*u*q+chi*q**2))

mu = 500
D = 50
chi = 1000
e0 = -9
rho = 530
ka = 500
kd =  0.0617

a_max = 9000
b_max = 1
c_max = 0.1

numN = 50
rng = 0.25

u = np.linspace(-rng,rng,numN)
q = np.linspace(-rng,rng,numN)
du = u[1]-u[0]

u, q = np.meshgrid(u, q)
flat = np.column_stack((u.flatten(), q.flatten()))


u_num = u.shape[0]
q_num = q.shape[0]

mapper = cm.ScalarMappable(Normalize(np.amin(0), np.amax(8000), clip=False), cm.plasma)

##### Plot initial surface #####
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax_a = plt.axes([0.25, 0.1, 0.65, 0.03])
a_slider = Slider(ax=ax_a, label=r'$a$', valinit=0, valstep=1, valmin=0, valmax=a_max)
a_slider.label.set_fontsize(24)

ax_b = plt.axes([0.25, 0.05, 0.65, 0.03])
b_slider = Slider(ax=ax_b, label=r'$b$', valinit=0, valstep=0.01, valmin=0, valmax=b_max)
b_slider.label.set_fontsize(24)

ax_c = plt.axes([0.25, 0.95, 0.65, 0.03])
c_slider = Slider(ax=ax_c, label=r'$c$', valinit=0, valstep=0.001, valmin=0, valmax=c_max)
c_slider.label.set_fontsize(24)


def update(val):

    a_val = a_slider.val
    b_val = b_slider.val
    c_val = c_slider.val
    
    n = dbb(u, q, mu, D, chi, e0, rho) - omega(u, q, a_val, b_val, c_val,  mu, D, chi)
    n /= (1+np.sum(n)*(u[1]-u[0])**4)
    nlt0 = np.where(n<-0.01)[0].shape[0]
    
    
    ex = dbb(u, q, mu, D, chi, e0, rho) 
    om = omega(u, q, a_val, b_val, c_val, mu, D, chi)
    
    intar = np.sum((ex-om)*du**2)
    print("Instances of negatives:  " + str(nlt0) + ". Integrated area:  " + str(intar) + ". Unbound fraction: " + str(1/(1+intar)))
    
    ax.clear() # Clear axes
    #ax.plot_surface(u, q, n, facecolors=mapper.to_rgba(n), lw=0, antialiased=False)
    #ax.scatter(u.flatten(), q.flatten(), ex.flatten(), c='crimson')
    #ax.scatter(u.flatten(), q.flatten(), om.flatten(), c='royalblue')
    ax.scatter(u.flatten(), q.flatten(), ex.flatten()-om.flatten(), c='k', alpha=0.5)
    ax.scatter([0], [0], [np.amax(n)], c='none')
    ax.set_xlabel('u', size=18)
    ax.set_ylabel('q', size=18)
    ax.set_zlabel('n', size=18)


a_slider.on_changed(update)
b_slider.on_changed(update)
c_slider.on_changed(update)

plt.show()