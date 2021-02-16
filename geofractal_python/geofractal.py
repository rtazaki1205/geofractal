import numpy as np
import math 
from scipy.special import gammainc, gammaincc

def geofractal(PN,df,k0,cormodel):
    
    """ 
    --------------------------------------------------------------------------------

       Python version of GEOFRACTAL version 1.0

    --------------------------------------------------------------------------------

     GEOFRACTAL computes geometrical cross sections of randomly orientated
     fractal dust aggregates based on a statistical distribution model of monomer
     particles developed in Tazaki (2021).

    --------------------------------------------------------------------------------
     INPUT PARAMETERS
    --------------------------------------------------------------------------------
     
     Df            : fractal dimension
     PN            : number of monomers
     k0            : fractal prefactor
     cormodel      : The cut-off model of the correlation function
                     cormodel = "EXPNL" : The exponential cut-off function
                     cormodel = "GAUSS" : The Gaussian cut-off function
                     cormodel = "FLDIM" : The fractal dimension cut-off function

    --------------------------------------------------------------------------------
     OUTPUT PARAMETER
    --------------------------------------------------------------------------------

     G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2.
          Therefore, "G" is non-dimensional.

    -------------------------------------------------------------------------------- 
    """

    #
    # safety checks
    #
    if cormodel not in 'EXPNL-GAUSS-FLDIM':
        print (' error: incorrect cormodel ')
        print (' stop ')
        exit()

    if PN < 0.9999:
        print (' error: number of monomer is less than 1.')
        print (' stop ')
        exit()

    if df < 0.9999 or df > 3.0001:
        print (' error: fractal dimension is out of its range.')
        print (' stop ')
        exit()

    #
    # Threshold number of monomers Nth
    #
    PNTH=min(11.0*df-8.5,8.0)

    #
    # calculation based on the analytical formula
    #
    if PN < PNTH:
        G = minato(PN)
    else:
        sigth = overlap(PNTH,k0,df,cormodel)
        Gth = minato(PNTH)
        A = (1.0+(PNTH-1.0)*sigth)*Gth
        sig = overlap(PN,k0,df,cormodel)
        G = A / (1.0+(PN-1.0)*sig)
    
    #
    # return the cross section
    #
    return G


def minato(PN):
    """ 
    Geometric cross sections of small non-fractal cluster    
    Equation (A.1) in Minato et al. (2006)
    """ 

    G = 12.5 * PN ** (-0.315) * np.exp(-2.53/PN**(0.0920))
    return G

def overlap(PN,k0,df,cormodel):
    """ 
    Calculate the overlapping efficiency
    """ 
  
    m = 999      # number of grids in numerical integration for a<0
    xmax = 25.0  # upper bound of numerical integration.
    
    if cormodel == 'EXPNL':
        a = df-2.0
        xmin = 2.0 * np.sqrt(0.5*df*(df+1.0))*(k0/PN)**(1.0/df)
        fac = xmin * xmin / 16.0 / math.gamma(df)
    elif cormodel == 'GAUSS':
        a = 0.5*(df-2.0)
        xmin = df*(k0/PN)**(2.0/df)
        fac = xmin / 16.0 / math.gamma(0.5*df)
    elif cormodel == 'FLDIM':
        a = (df-2.0)/df
        xmin = 2.0 ** (df-1.0) * k0 / PN
        fac = xmin ** (2.0/df)/16.0

    if a > 0:
        intg = math.gamma(a)*gammaincc(a,xmin)
    else:
        x = np.exp(np.linspace(math.log(xmin),math.log(xmax),m+1))
        dlnx = math.log(xmax/xmin)/m
        f = x**a*np.exp(-x)
        intg = 0.0 
        for i in range(m):
            intg += 0.5*(f[i]+f[i+1])
        intg = intg * dlnx

    sig = fac*intg

    return sig 
