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

def integrand(om,u,q,du,dq,mu,D,chi,rho,e0,ka,kd,a,b,c):
    pkbt = rho*1.381E-23*300
    ex = np.exp(-e0-0.5*(mu*u**2+2*D*u*q+chi*q**2)/pkbt)-a*np.exp(-b*((u-c)**2+q**2))
    
    #arr = np.exp(-0.5*3E-4*(u*40E-9)**2/(300*1.381E-23))
    #ka_ar = ka*arr
   
    al1 = (1+np.sum(ex*du*dq))
    mu_u = mu*u
    d_q = D*q
    nu_bar = kd*(ex)/(1+ex)
    return (4/al1)*ex*(mu_u**2 + 2*mu_u*d_q + d_q**2)*nu_bar/(nu_bar**2+(2*np.pi*om)**2)


n = 10
mapper = cm.ScalarMappable(Normalize(0,8100), cm.plasma)
a_ar = []
plateau = []
unbd = []
al1_ar = []

plt.figure()
 
for i in range(n):
    fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Report1_0', '1d_uq_om_') + str(i)
    psd = np.loadtxt(fold+"/PSD.dat") # t, stress, orientation
    
    pop_distrib = np.loadtxt(fold+"/N_Final.dat")
    unbd_whr = np.where(pop_distrib[:,0] == -2)[0]
    unbd.append(np.mean(pop_distrib[unbd_whr,1])/1.224490e-02)
    
    
    f = open(fold + '/parameters.dat')

    for l in f:
        
        if l[:2] == 'mu':
            mu = float(l[4:])
        elif l[:2] == 'D:':
            D_qu = float(l[3:])
        elif l[:3] == 'chi':
            chi = float(l[5:])
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
    
    a_ar.append(a)
    
    colour = mapper.to_rgba(a)
    # plt.errorbar(psd[:,0], psd[:,1], yerr=psd[:,2], markeredgecolor=colour, ecolor=colour,
    #               markerfacecolor='none', markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1)
    plt.scatter(psd[:,0], 2*np.pi*psd[:,1]*psd[:,0], color=colour, fc='none', s=70, marker='v')
    
    om = np.logspace(-6,np.log10(50),100)
    u = np.linspace(-0.25,0.25,1000)
    q = np.linspace(-0.25,0.25,1000)

    u_m, q_m = np.meshgrid(u, q)
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    dq = (np.amax(q)-np.amin(q))/(q.shape[0]-1)
    sig = []
    for o in om:
        sig.append(np.sum(integrand(o, u_m, q_m, du, dq, mu, D_qu, chi, rho, e0, ka, kd, a, b, c)*du*dq))
    plt.plot(om, np.array(sig)*2*np.pi*om, c=colour, alpha=0.5, lw=3)
    
    plateau.append(sig[0])
    
    pkbt = rho*1.381E-23*300
    ex = np.exp(-e0-0.5*mu*u**2/pkbt)-a*np.exp(-b*(u-c)**2)
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    al1 = (1+np.sum(ex*du))
    al1_ar.append(al1)
    
    #plt.plot(u, a*np.exp(-b*(u-c)**2), color=colour, lw=2)
    
plt.xlim([1E-5, 10])
plt.ylim([0.01,5E5])   
    
cb = plt.colorbar(mapper, label='c')
cb.set_label(label=r'$c$', size='xx-large')
cb.ax.tick_params(labelsize='large')
    
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$f/\textrm{Hz}$", size=24)
plt.ylabel(r"$|\sigma(f)|^2 / \textrm{Pa}^2\textrm{s}$", size=24)
plt.ylabel(r"$\alpha^{\prime\prime},\, \omega S(\omega) / 2k_B T (\textrm{Pa}^2)$", size=24)
plt.tick_params(which='major', direction='in', width=3, labelsize=20, length=6, pad=5, top=True, right=True)
plt.tick_params(which='minor', direction='in', width=1.5, length=3, top=True, right=True)
    
    
    
    
    
    
    
    
    
    
    