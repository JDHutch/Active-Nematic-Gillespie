import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from scipy.optimize import curve_fit as cf

## Latex Fonts ##
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('legend', fancybox=False, edgecolor='k', shadow=False)
plt.rc('axes', linewidth=3)

def integrand(om,u,q,du,dq,mu,D,chi,rho,e0,ka,kd,a,b,c):
    pkbt = rho*1.381E-23*300
    #ex = np.exp(-e0-0.5*(mu*u**2+2*D*u*q+chi*q**2)/pkbt)-a*np.exp(-b*((u-c)**2+q**2))
    ex = np.exp(-e0-0.5*(mu*u**2+2*D*u*q+chi*q**2)/pkbt)-a*np.exp(-0.5*b*(mu*(u-c)**2+2*D*q*u+chi*q**2))
    
    #arr = np.exp(-0.5*3E-4*(u*40E-9)**2/(300*1.381E-23))
    #ka_ar = ka*arr
   
    al1 = (1+np.sum(ex*du*dq))
    mu_u = mu*u
    d_q = D*q
    nu_bar = kd*(ex)/(1+ex)
    return (4/al1)*ex*(mu_u**2 + 2*mu_u*d_q + d_q**2)*nu_bar/((kd+ kd/al1)**2+(2*np.pi*om)**2)


### Fitting ####
def integrand_fit(om,u,q,du,dq,A):
    mu = 500
    D = 50
    chi = 1000
    kd = 0.062
    a = 8100
    b = 0.46
    c = 0
    e0 = -9
    
    
    pkbt = 530E18*1.381E-23*300
    ex = np.exp(-e0-0.5*(mu*u**2+2*D*u*q+chi*q**2)/pkbt)-a*np.exp(-0.5*b*(mu*(u-c)**2+2*D*q*u+chi*q**2))
    
    ka_ar = kd*ex
    

    al1 = (1+np.sum(ex*du*dq))
    mu_u = mu*u
    d_q = D*q
    ka_ar /= al1

    nu_bar = kd*ka_ar/(kd+ka_ar)
    return (4/al1)*ex*(mu_u**2 + 2*mu_u*d_q + d_q**2)*nu_bar/((kd/A)**2+(2*np.pi*om)**2)

def igft(om, A):
    sig = []
    u = np.linspace(-0.25,0.25,1000)
    q = np.linspace(-0.25,0.25,1000)
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    dq = (np.amax(q)-np.amin(q))/(q.shape[0]-1)
    u, q = np.meshgrid(u, q)
    for o in om:
        sig.append(np.sum(integrand_fit(o,u, q, du, dq, A)*du*dq))
    return np.array(sig)

#################


n = 10
mapper = cm.ScalarMappable(Normalize(0,8100), cm.viridis)
a_ar = []
al1_ar = []
fits = []


plt.figure()
 
for i in range(n):
    fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Report1_1', '1d_uq_om_') + str(i)
    psd = np.loadtxt(fold+"/PSD.dat") # t, stress, orientation 
    
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
    plt.errorbar(psd[:,0], psd[:,1], yerr=psd[:,2], markeredgecolor=colour, ecolor=colour,
                  markerfacecolor='none', markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1)
    #plt.scatter(psd[:,0], 2*np.pi*psd[:,1]*psd[:,0], color=colour, fc='none', s=70, marker='v')

    ft = cf(igft, psd[:,0], psd[:,1])
    fits.append(ft[0])


    plt.plot(psd[:,0], igft(psd[:,0], ft[0]), c=colour, alpha=0.5, lw=3)
    
    
    u = np.linspace(-0.25,0.25,1000)
    q = np.linspace(-0.25,0.25,1000)
    u_m, q_m = np.meshgrid(u, q)
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    dq = (np.amax(q)-np.amin(q))/(q.shape[0]-1)
    pkbt = rho*1.381E-23*300
    ex = np.exp(-e0-0.5*(mu*u_m**2+2*D_qu*q_m*u_m+chi*q_m**2)/pkbt)-a*np.exp(-0.5*b*(mu*(u_m-c)**2+2*D_qu*q_m*u_m+chi*q_m**2))
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    al1 = (1+np.sum(ex*du**2))
    al1_ar.append(al1)
    

# plt.xlim([1E-5, 1])
# plt.ylim([1,7.5E8])   
    
cb = plt.colorbar(mapper, label='a')
cb.set_label(label=r'$a$', size='xx-large')
cb.ax.tick_params(labelsize='large')
    
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$f/\textrm{Hz}$", size=24)
plt.ylabel(r"$|\sigma(f)|^2 / \textrm{Pa}^2\textrm{s}$", size=24)
#plt.ylabel(r"$\alpha^{\prime\prime},\, \omega S(\omega) / 2k_B T (\textrm{Pa}^2)$", size=24)
plt.tick_params(which='major', direction='in', width=3, labelsize=20, length=6, pad=5, top=True, right=True)
plt.tick_params(which='minor', direction='in', width=1.5, length=3, top=True, right=True)
    
    
    
    
    
    
    
    
    
    
    