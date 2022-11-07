import matplotlib
font = {'family' : 'Tahoma',
        'weight' : 'bold',
        'size'   : 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import forward
from math import exp
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import loadmat
from random import randint, gauss, random as rand
import time





rho_new = 300.
eta = .9*10**7*86400.
depth0 = 6.120-RHdata[0]
#acc0 =.1*np.random.rand(200)
acc0 = 0.3*np.ones(200)
#acc0=np.diff(RH_data['RHm200'],axis=0)
#acc0=np.append([0.],-1.0*acc0)
#acc0[acc0<0]=0.0
D0=depth0/40.*np.ones(40, dtype=float)
rho0=np.arange(300,400,100./len(D0))
l=len(acc0)+len(rho0)+2

outputs, misfits = forward.monte_carlo_run(rho0,acc0,rho_new,eta,D0,RH_data['RHm200'],n_iterates=20000,regularization=0.0)#1000.0)



best = np.argmin(misfits[:, 0])
rho0best=outputs[:,best][:len(rho0)]
accbest=outputs[:,best][len(rho0):][:len(acc0)]
rho_newbest=outputs[:,best][len(rho0):][len(acc0):][0]
etabest=outputs[:,best][len(rho0):][len(acc0):][1]
print(rho0best)
print(etabest)
print(rho_newbest)

doy200 = np.arange(0,200,1)
depth = np.cumsum(D0)

fig, ax = plt.subplots(figsize=(12, 7))
plt.plot(np.arange(misfits.shape[0]), misfits[:, 0])
ax.set_xlabel('Iterate')
ax.set_ylabel('Cost')
plt.savefig('model_misfit.png')
plt.close(fig)


fig2, ax2 = plt.subplots(figsize=(12, 7))
RHbest = forward.model(D0, accbest,rho0best,200*86400,rho_newbest,etabest,n_tsteps=200)
ax2.plot(doy200, RH_data['RHm200'])
ax2.plot(doy200,RHbest , color='k')
ax2.set_xlabel('day of year')
ax2.set_ylabel('receiver height.(m)')
plt.savefig('best_solution.png')
plt.close(fig)


# Plot it---it is linear as we expect
plt.figure(figsize=(12, 7))
plt.title('Vertical density profile')
plt.plot(rho0best, depth)
plt.xlabel('rho.(kg/m^3)')
plt.ylabel( 'depth.(m)')
plt.gca().invert_yaxis()
plt.savefig('density.png')
plt.close(fig)

plt.figure(figsize=(12,7))
plt.plot(doy200,accbest)
plt.title('Accumulation rate.(m/d)')
plt.xlabel('day of year')
plt.ylabel('accumulation.(m)')
plt.savefig('accumulation.png')
plt.close(fig)
