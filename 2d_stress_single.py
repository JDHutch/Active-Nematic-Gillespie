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

def integrand(om,ua,ub,du,mu,rho,kbT,e0,ka,kd,a,b,c):
    pkbt = rho*kbT*1E18*1.381E-23*300
    ex = np.exp(-e0/kbT-0.5*mu*(ua**2+ub**2)/pkbt)-a*np.exp(-b*((ua-c)**2+(ub-c)**2))
    #omega = a*np.exp(-b*(u-c)**2)
    #arr = np.exp(-0.5*3E-4*(u*40E-9)**2/(300*1.381E-23))
    #ka_ar = ka*arr
    al1 = (1+np.sum(ex*du*du))
    return 2*(mu**2/al1)*u_a**2*kd*ex/((kd)**2+(2*np.pi*om)**2)  #2*(mu**2/al1)*u**2/((ka_ar/ex)**2 + (2*np.pi*om)**2)

def ashot(om, u, mu, rho, kbT, e0, ka, a, b, c):
    pkbt = rho*kbT*1E18*1.381E-23*300
    ex = np.exp(-e0/kbT-0.5*mu*u**2/pkbt)-a*np.exp(-b*(u-c)**2)
    arr = np.exp(-0.5*3E-4*(u*40E-9)**2/(300*1.381E-23))
    ka_ar = ka*arr
    du = (np.amax(u)-np.amin(u))/(u.shape[0]-1)
    al1 = (1+np.sum(ex*du))
    return 4*np.pi*mu*pkbt*ka_ar/(al1)/((ka_ar/ex)**2 + (2*np.pi*om)**2)

def int_mul(y,a):
    return a*y

### Load file ###
fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Tau_test', '2d_tau')
stress = np.loadtxt(fold+"/stress_0.dat") # t, stress, orientation
#stress[:,1] -= np.mean(stress[:,1])

# Sampling rate
dt = 0.005 # np.mean(stress[1:,0]-stress[:-1,0])

## Constants ##

X_T = 1000000

mu = 500
rho = 530
kbT = 1
e0 = -9
ka = 500
kd =  0.0617

a = 0
b = 150
c = 0

tunq, unqwhr = np.unique(stress[:,0], return_index=True)
stress = stress[unqwhr,:]


interp = interp1d(stress[:,0], stress[:,1], kind='cubic')
interp2 = interp1d(stress[:,0], stress[:,2], kind='cubic')
newt = np.arange(stress[0,0],stress[-1,0],dt)
newstress = interp(newt)
newstress2 = interp2(newt)

f = np.fft.fftfreq(newt.shape[0], dt) # Freqs
ps = np.fft.fft(newstress) # Power spectrum
ps2 = np.fft.fft(newstress2) # Power spectrum

n_oneside = newt.shape[0]//2
ps = abs(ps[:n_oneside]*dt)**2/(tunq[-1]-tunq[0]) ## Normalise for sample length and sampling rate
ps2 = abs(ps2[:n_oneside]*dt)**2/(tunq[-1]-tunq[0])
f = f[:n_oneside]


pMean, fLog = bs(f, ps, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='mean')[:2]
pStd = bs(f, ps, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='std')[0]
pStd /= np.sqrt(bs(f, ps, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='count')[0])
pMean2 = bs(f, ps2, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='mean')[0]
pStd2 = bs(f, ps2, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='std')[0]
pStd2 /= np.sqrt(bs(f, ps2, bins = np.logspace(np.log10(f[1]), np.log10(f[-1]),50), statistic='count')[0])

pNan2 = np.isnan(pMean2)
pMean2 = pMean2[~pNan2]
pStd2 = pStd2[~pNan2]
fLog2 = 0.5*(fLog[1:]+fLog[:-1])[~pNan2]
pNan = np.isnan(pMean)
pMean = pMean[~pNan]
pStd = pStd[~pNan]
fLog = 0.5*(fLog[1:]+fLog[:-1])[~pNan]

colour = 'crimson'

#pStd[5] *= 0.75

#plt.figure()
plt.errorbar(fLog, pMean, yerr=pStd, markeredgecolor=colour, ecolor=colour,
             markerfacecolor='none', markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1, label=r"$\textrm{Tau-Leaping: }\sigma_{xx}=-\sigma_{yy}$")
plt.errorbar(fLog, pMean2, yerr=pStd, markeredgecolor=colour, ecolor=colour,
             markerfacecolor=colour, markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1, label=r"$\textrm{Tau-Leaping: }\sigma_{xy}=\sigma_{yx}$")
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$f/\textrm{Hz}$", size=24)
plt.ylabel(r"$|\sigma(f)|^2 / \textrm{Pa}^2\textrm{s}$", size=24)
plt.tick_params(which='major', direction='in', width=3, labelsize=20, length=6, pad=5, top=True, right=True)
plt.tick_params(which='minor', direction='in', width=1.5, length=3, top=True, right=True)

om = np.logspace(-6,np.log10(10),100)
ua = np.linspace(-1,1,100)
ub = np.linspace(-1,1,100)
du = (np.amax(ua)-np.amin(ua))/(ua.shape[0]-1)
u_a, u_b = np.meshgrid(ua,ub)
sig = []
for o in om:
    sig.append(np.sum(integrand(o, u_a, u_b, du, mu, rho, kbT, e0, ka, kd, a, b, c)*du**2))
plt.plot(om, np.array(sig), c='k', alpha=0.5, lw=2)

om = np.logspace(-6,np.log10(1),100)
u = np.linspace(-1,1,100)
sig_eq = []
for o in om:
    sig_eq.append(np.trapz(ashot(o, u, mu, rho, kbT, e0, ka,a,b,c),u))
#plt.plot(om, np.array(sig_eq), c='k', alpha=0.5, lw=2)

plt.xlim([7E-6, 2])
plt.ylim([6,1E5])