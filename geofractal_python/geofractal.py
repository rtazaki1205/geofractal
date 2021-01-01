import numpy as np
import math 
from scipy.special import gammainc, gammaincc

def geofractal(PN,df,k0):
    if df > 1.5:
        PNTH = 8.0
    else:
        PNTH = 11.0 * (df-1.0) + 2.5

    if PN < PNTH:
        G = minato(PN)
        return G
    else:
        sigth = overlap(PNTH,k0,df)
        Gth = minato(PNTH)
        A = (1.0+(PNTH-1.0)*sigth)*Gth
        sig = overlap(PN,k0,df)
        G = A / (1.0+(PN-1.0)*sig)
        return G

def minato(PN):
    G = 12.5 * PN ** (-0.315) * np.exp(-2.53/PN**(0.0920))
    return G

def overlap(PN,k0,df):
    a = (df-2.0)/df
    eta = 2.0 ** (df-1.0) * k0 / PN

    if a > 0:
        intg = math.gamma(a)*gammaincc(a,eta)
    else:
        m = 1000
        xmin = eta
        xmax = 25.0
        x = np.exp(np.linspace(math.log(xmin),math.log(xmax),m+1))
        dlnx = math.log(xmax/xmin)/(m+1)
        f = x**a*np.exp(-x)
        intg = 0.0 
        for i in range(m):
            intg += 0.5*(f[i]+f[i+1])
        intg = intg * dlnx

    sig = eta**(2.0/df)*intg/16.0
    return sig 
