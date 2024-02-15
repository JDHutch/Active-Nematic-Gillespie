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

def dbb(ua, ub, qa, qb, mu, D, chi, e0, rho):
    xx = mu*ua**2+2*D*ua*qa+chi*qa**2
    xy = mu*ub**2+2*D*ub*qb+chi*qb**2
    return np.exp(-e0-(xx+xy)/(rho*1E18*300*1.381E-23))

def omega(ua, ub, qa, qb, a, b, c):
    return a*np.exp(-b*((ua-c)**2+ub**2+qa**2+qb**2))

mu = 500
D = 50
chi =  1000
e0 = -9
rho = 530
ka = 500
kd =  0.0617

a = 0
b = 1500
c = 0.01

numN = 25
rng = 0.4

ua = np.linspace(-rng,rng,numN)
ub = np.linspace(-rng,rng,numN)
qa = np.linspace(-rng,rng,numN)
qb = np.linspace(-rng,rng,numN)

u_a, u_b, q_a, q_b = np.meshgrid(ua, ub, qa, qb)
flat = np.column_stack((u_a.flatten(), u_b.flatten(), q_a.flatten(), q_b.flatten()))

qu_plots = [ua, ub, qa, qb]

n = dbb(flat[:,0], flat[:,1], flat[:,2], flat[:,3], mu, D, chi, e0, rho) - omega(flat[:,0], flat[:,1], flat[:,2], flat[:,3], a, b, c)
n /= (1+np.sum(n)*(ua[1]-ua[0])**4)
nlt0 = np.where(n<-0.01)[0].shape[0]
print("Instances of negatives: " + str(nlt0))

ua_num = ua.shape[0]
ub_num = ua.shape[0]
qa_num = qa.shape[0]
qb_num = qb.shape[0]

mapper = cm.ScalarMappable(Normalize(np.amin(n), np.amax(n), clip=False), cm.plasma)

##### Plot initial surface #####
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax_ua = plt.axes([0.25, 0.1, 0.65, 0.03])
ua_slider = Slider(ax=ax_ua, label=r'$u_{a}\,\textrm{Slice}$', valinit=0, valstep=1, valmin=0, valmax=ua_num-1)
ua_slider.label.set_fontsize(24)

ax_ub = plt.axes([0.25, 0.05, 0.65, 0.03])
ub_slider = Slider(ax=ax_ub, label=r'$u_{b}\,\textrm{Slice}$', valinit=0, valstep=1, valmin=0, valmax=ub_num-1)
ub_slider.label.set_fontsize(24)

ax_qa = plt.axes([0.25, 0.95, 0.65, 0.03])
qa_slider = Slider(ax=ax_qa, label=r'$q_{a}\,\textrm{Slice}$', valinit=0, valstep=1, valmin=0, valmax=qa_num-1)
qa_slider.label.set_fontsize(24)

ax_qb = plt.axes([0.25, 0.9, 0.65, 0.03])
qb_slider = Slider(ax=ax_qb, label=r'$q_{b}\,\textrm{Slice}$', valinit=0, valstep=1, valmin=0, valmax=qb_num-1)
qb_slider.label.set_fontsize(24)

ax_r1 = plt.axes([0.075, 0.5, 0.1, 0.3])
r1 = RadioButtons(ax_r1, (r'$u_{a}$', r'$u_{b}$', r'$q_{a}$', r'$q_{b}$'), active=0, activecolor='k')
for j in r1.labels:
    j.set_fontsize(24)

ax_r2 = plt.axes([0.075, 0.2, 0.1, 0.3])
r2 = RadioButtons(ax_r2, (r'$u_{a}$', r'$u_{b}$', r'$q_{a}$', r'$q_{b}$'), active=0, activecolor='k')
for j in r2.labels:
    j.set_fontsize(24)


radio_dict = {r'$u_{a}$': 0, r'$u_{b}$': 1, r'$q_{a}$': 2, r'$q_{b}$': 3}

x_ax = 0
y_ax = 1

def update(val):

    x_ax = radio_dict[r1.value_selected]
    y_ax = radio_dict[r2.value_selected]

    cnst = list(set([0,1,2,3])-set([x_ax,y_ax]))

   # ua_i = ua_slider.val
    #ub_i = ub_slider.val
    #qa_i = qa_slider.val
    #qb_i = qb_slider.val
    
    slds = [ua_slider.val, ub_slider.val, qa_slider.val, qb_slider.val]
   
    whr = np.where((flat[:,cnst[0]] == ua[slds[cnst[0]]])*(flat[:,cnst[1]] == ua[slds[cnst[1]]]))[0]
    
    n_plot = n[whr].reshape(numN,numN)
    
    x = flat[:,x_ax][whr].reshape(numN,numN)
    y = flat[:,y_ax][whr].reshape(numN,numN)
    
    ax.clear() # Clear axes
    ax.plot_surface(x, y, n_plot, facecolors=mapper.to_rgba(n_plot), lw=0, antialiased=False)
    #ax.scatter(x.flatten(), y.flatten(), n_plot.flatten(), c=mapper.to_rgba(N.flatten()))
    ax.scatter([0], [0], [np.amax(n)], c='none')
    ax.set_xlabel(r1.value_selected, size=18)
    ax.set_ylabel(r2.value_selected, size=18)
    ax.set_zlabel('n', size=18)


ua_slider.on_changed(update)
ub_slider.on_changed(update)
qa_slider.on_changed(update)
qb_slider.on_changed(update)

r1.on_clicked(update)
r2.on_clicked(update)

plt.show()