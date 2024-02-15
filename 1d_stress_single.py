import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic as bs
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit as cf

## Latex Fonts ##
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('legend', fancybox=False, edgecolor='k', shadow=False)
plt.rc('axes', linewidth=3)

def integrand(om,u,mu,rho,e0,ka,kd,a,b,c):
    pkbt = rho*1.381E-23*300
    ex = np.exp(-e0-0.5*mu*u**2/pkbt)-a*np.exp(-b*(u-c)**2)
    #omega = a*np.exp(-b*(u-c)**2)
    #arr = np.exp(-0.5*3E-4*(u*40E-9)**2/(300*1.381E-23))
    #ka_ar = ka*arr
    #kd = ka/ex
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    al1 = (1+np.sum(ex*du))
    nu_bar = kd*(ex)/(1+ex)
    return 4*(mu**2/al1)*u**2*nu_bar*ex/((kd)**2+(2*np.pi*om)**2)

def int_mul(y,a):
    return a*y

### Load file ###
#fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Report7', '1d_u_om_9') 
fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Tau_test', 'normal_tau') 
stress = np.loadtxt(fold+"/stress_0.dat") # t, stress, orientation
#stress[:,1] -= np.mean(stress[:,1])

## Constants ##

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



## Sampling rate ##
dt = np.mean(stress[1:,0]-stress[:-1,0])

tunq, unqwhr = np.unique(stress[:,0], return_index=True)
stress = stress[unqwhr,:]


interp = interp1d(stress[:,0], stress[:,1], kind='cubic')
newt = np.arange(stress[0,0],stress[-1,0],dt)
newstress = interp(newt)


f = np.fft.fftfreq(newt.shape[0], dt) # Freqs
ps = np.fft.fft(newstress) # Power spectrum

n_oneside = newt.shape[0]//2
ps = 2*abs(ps[:n_oneside]*dt)**2/(tunq[-1]-tunq[0]) ## Normalise for sample length and sampling rate x2 for positive freq only
f = f[:n_oneside]


pMean, fLog = bs(f, ps, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='mean')[:2]
pStd = bs(f, ps, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='std')[0]
pStd /= np.sqrt(bs(f, ps, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='count')[0])
pNan = np.isnan(pMean)
pMean = pMean[~pNan]
pStd = pStd[~pNan]
fLog = 0.5*(fLog[1:]+fLog[:-1])[~pNan]

colour = 'crimson'

plt.figure()
plt.errorbar(fLog, pMean, yerr=pStd, markeredgecolor=colour, ecolor=colour,
             markerfacecolor='none', markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1, label="Tau-leaping")
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$f/\textrm{Hz}$", size=24)
plt.ylabel(r"$|\sigma(f)|^2 / \textrm{Pa}^2\textrm{s}$", size=24)
plt.tick_params(which='major', direction='in', width=3, labelsize=20, length=6, pad=5, top=True, right=True)
plt.tick_params(which='minor', direction='in', width=1.5, length=3, top=True, right=True)

om = np.logspace(-6,np.log10(50),100)
u = np.linspace(-0.25,0.25,100)
sig = []
for o in om:
    sig.append(np.trapz(integrand(o, u, mu, rho, e0, ka, kd, a, b, c),u))
plt.plot(om, np.array(sig), c='crimson', alpha=0.5, lw=2)

sig_ft = []
for o in fLog:
    sig_ft.append(np.trapz(integrand(o, u, mu, rho, e0, ka, kd, a, b, c),u))
mult = cf(int_mul, np.array(sig_ft[:]), pMean[:])[0]
#plt.plot(om, np.array(sig)*mult[0], c='royalblue', lw=2)
print("Multiplier: " + str(mult[0]))


plt.xlim([0.5E-6, 20])
plt.ylim([0.02,1E6])

np.savetxt(fold + '/PSD.dat', np.column_stack((fLog, pMean, pStd)))
