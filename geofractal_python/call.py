import numpy as np
import math
from geofractal import *

df = 3.0
k0 = 0.5*(0.3-np.sqrt(3.0))*(df-1.0)+np.sqrt(3.0)

Nmin = 1.e0
Nmax = 1.e10
N    = 250
PN = np.exp(np.linspace(math.log(Nmin),math.log(Nmax),N))

G = np.zeros(N)
for i in range(N):
    G[i] = geofractal(PN[i],df,k0)

filename='gratio.out'
with open(filename,'w') as f:
    f.write('# df  = %13.6e \n'%df)
    f.write('# k0  = %13.6e \n'%k0)
    f.write('# %11s %13s\n'%('N','G/NpiR0^2'))
    for i in range(N):
        f.write('%13.6e %13.6e\n'%(PN[i],G[i]))
