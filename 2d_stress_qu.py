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

def integrand(om,u_a,u_b,q_a,q_b,du,dq,mu,D,chi,rho,kbT,e0,ka,kd,a,b,c):
    pkbt = rho*kbT*1E18*1.381E-23*300
    ex = np.exp(-e0/kbT-(mu*(u_a**2+u_b**2)+2*D*(u_a*q_a+u_b*q_b)+chi*(q_a**2+q_b**2))/pkbt)
    ex -= a*np.exp(-b*((u_a-c)**2 + u_b**2 + q_a**2 + q_b**2))
    #arr = np.exp(-0.5*3E-4*(u*40E-9)**2/(300*1.381E-23))
    #ka_ar = ka*arr
    
   
    al1 = (1+np.sum(ex*du*dq*du*dq))
    mu_u_a = mu*u_a
    d_q_a = D*q_a
    nu_bar = kd*(ex)/(1+ex)
    return (4/al1)*ex*(mu_u_a**2 + 2*mu_u_a*d_q_a + d_q_a**2)*nu_bar/(nu_bar**2+(2*np.pi*om)**2)

def int_mul(y,a):
    return a*y

### Load file ###
fold = os.path.join(os.path.expanduser('~'), 'Documents', 'C', 'Gillespie_C', 'Tau_test', '2d_qu_tau_7')
stress = np.loadtxt(fold+"/stress_0.dat") # t, stress, orientation
#stress[:,1] -= np.mean(stress[:,1])
#stress = stress[200000:,:]
#stress[:,0] -= stress[0,0]

# Sampling rate
dt = 0.001# # np.mean(stress[1:,0]-stress[:-1,0])

## Constants ##

X_T = 1000000

mu = 1000
D_qu = 50
chi = 1000
rho = 530
kbT = 1
e0 = -9
ka = 500
kd =  0.0617

a = 0
b = 1500
c = 0.01

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
ps = 2*abs(ps[:n_oneside]*dt)**2/(tunq[-1]-tunq[0]) ## Normalise for sample length and sampling rate
ps2 = 2*abs(ps2[:n_oneside]*dt)**2/(tunq[-1]-tunq[0])
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

colour = 'royalblue'

#pStd[5] *= 0.75

plt.figure()
plt.errorbar(fLog, pMean, yerr=pStd, markeredgecolor=colour, ecolor=colour,
             markerfacecolor='none', markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1, label=r"$\textrm{Gillespie: }\sigma_{xx}=-\sigma_{yy}$")
plt.errorbar(fLog, pMean2, yerr=pStd, markeredgecolor=colour, ecolor=colour,
             markerfacecolor=colour, markersize=8, fmt='v', capsize=5, capthick=1, elinewidth=1, label=r"$\textrm{Gillespie: }\sigma_{xy}=\sigma_{yx}$")
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$f/\textrm{Hz}$", size=24)
plt.ylabel(r"$|\sigma(f)|^2 / \textrm{Pa}^2\textrm{s}$", size=24)
plt.tick_params(which='major', direction='in', width=3, labelsize=20, length=6, pad=5, top=True, right=True)
plt.tick_params(which='minor', direction='in', width=1.5, length=3, top=True, right=True)

om = np.logspace(-6,np.log10(10),100)
ua = np.linspace(-1,1,50)
qa = np.linspace(-1,1,50)
ub = np.linspace(-1,1,50)
qb = np.linspace(-1,1,50)
u_ma, u_mb, q_ma, q_mb = np.meshgrid(ua, ub, qa, qb)
du = (np.amax(ua)-np.amin(ua))/(ua.shape[0]-1)
dq = (np.amax(qa)-np.amin(qa))/(qa.shape[0]-1)
sig_a = []
sig_b = []
for o in om:
    sig_a.append(np.sum(integrand(o, u_ma, u_mb, q_ma, q_mb, du, dq, mu, D_qu, chi, rho, kbT, e0, ka, kd, a, b, c)*du*dq*du*dq))
#     sig_b.append(np.sum(integrand(o, u_mb, q_mb, du, dq, mu, D_qu, chi, rho, kbT, e0, ka, kd, a, b, c)*du*dq))
plt.plot(om, np.array(sig_a), c='k', alpha=0.5, lw=2)

sig_ft = []
for o in fLog:
    sig_ft.append(np.sum(integrand(o, u_ma, u_mb, q_ma, q_mb, du, dq, mu, D_qu, chi, rho, kbT, e0, ka, kd, a, b, c)*du*dq*du*dq))
mult = cf(int_mul, np.array(sig_ft[15:]), pMean[15:])[0]
plt.plot(om, np.array(sig_a)*mult[0], c='royalblue', lw=2)
print("Multiplier: " + str(mult[0]))



plt.xlim([7E-6, 2])
#plt.ylim([0.001,1E10])