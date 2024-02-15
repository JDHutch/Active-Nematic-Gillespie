import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

## Latex Fonts ##
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('legend', fancybox=False, edgecolor='k', shadow=False)
plt.rc('axes', linewidth=3)

def integrand(om,u,mu,rho,e0,ka,kd,a,b,c):
    pkbt = rho*1.381E-23*300
    ex = np.exp(-e0-0.5*mu*u**2/pkbt)-a*np.exp(-b*(u-c)**2)
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    al1 = (1+np.sum(ex*du))
    nu_bar = kd*(ex)/(1+ex)
    return 4*(mu**2/al1)*u**2*nu_bar*ex/((kd)**2+(2*np.pi*om)**2)


n = 10
mapper = cm.ScalarMappable(Normalize(0,n-1), cm.hsv)

i = 0
 
fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Report', '1d_u_om_') + str(i)
psd = np.loadtxt(fold+"/PSD.dat") # t, stress, orientation

f = open(fold + '/parameters.dat')

for l in f:
    
    if l[:2] == 'mu':
        mu = float(l[4:])
    elif l[:2] == 'e0':
        e0 = float(l[4:])
    elif l[:3] == 'rho':
        rho = float(l[5:])
    elif l[:5] == 'k_a^0':
        ka =  float(l[7:])
    elif l[:3] == 'k_d':
        kd = float(l[5:])
    elif l[:9] == 'Prefactor':
        a = float(l[11:])
    elif l[:8] == 'Exponent':
        b = float(l[10:])
    elif l[:5] == 'Shift':
        c = float(l[7:])
f.close()

fig, ax1 = plt.subplots()
left, bottom, width, height = [0.275, 0.275, 0.375, 0.375]
ax2 = fig.add_axes([left, bottom, width, height])


colour = mapper.to_rgba(i)
ax1.errorbar(psd[:,0], psd[:,1], yerr=psd[:,2], markeredgecolor=colour, ecolor=colour,
             markerfacecolor='none', markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1)

om = np.logspace(-6,np.log10(50),100)
u = np.linspace(-0.25,0.25,100)
sig = []
for o in om:
    sig.append(np.trapz(integrand(o, u, mu, rho, e0, ka, kd, a, b, c),u))
ax1.plot(om, np.array(sig), c=colour, alpha=0.5, lw=2)

u = np.linspace(-0.25,0.25,1000)
pkbt = rho*1.381E-23*300
ex = np.exp(-e0-0.5*mu*u**2/pkbt)-a*np.exp(-b*(u-c)**2)
du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
al1 = (1+np.sum(ex*du))

ax2.plot(u, ex/al1, color=colour, lw=1.5)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_xlabel(r"$f/\textrm{Hz}$", size=24)
ax1.set_ylabel(r"$|\sigma(f)|^2 / \textrm{Pa}^2\textrm{s}$", size=24)
ax1.tick_params(which='major', direction='in', width=3, labelsize=20, length=6, pad=5, top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1.5, length=3, top=True, right=True)

ax2.set_xlabel(r"$u$", size=16, labelpad = 0.1)
ax2.set_ylabel(r"$n(u)$", size=16, labelpad = 0.1)
ax2.tick_params(which='major', direction='in', width=3, labelsize=16, length=3, pad=3, top=True, right=True)
    
ax1.set_title("a = " + str(a) + ", b = " + str(b) + ", c = " + str(c), size=24)
    
fold2 = os.path.join(os.path.expanduser('~'), 'Documents', 'Reports', 'Images')
file_name = fold2+'/PSD_Example' + str(i) + '.pdf'
if(False):
    plt.savefig(file_name)
    
    
    
    
    
    
    