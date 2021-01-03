import numpy as np

def minato06_bpca(PN):

    if PN > 16:
        pca = 4.27*PN**(-0.315)*np.exp(-1.74/PN**0.243)
    else:
        pca = 12.5*PN**(-0.315)*np.exp(-2.53/PN**0.0920)

    return pca

def minato06_bcca(PN):
    if PN > 16:
        cca = 0.352 + 0.566*PN**(-0.138)
    else:
        cca = 12.5*PN**(-0.315)*np.exp(-2.53/PN**0.0920)

    return cca

def chain_exact(PN):

    G = 1.0/PN + (1.0-1.0/PN)*8.0/np.pi/3.0

    return G


def chararea(PN,df,k0):
    G = 5.0 * PN ** (-1.0+2.0/df)/3.0/k0**(2.0/df)
    G = min(G,1.0)
    return G

def okuzumi09(PN,df,k0):
    k0_bcca = 1.04
    df_bcca = 1.9
    rc  = np.sqrt(5.0/3.0) * (PN/k0) ** (1.0/df)
    rc_bcca = np.sqrt(5.0/3.0) * (PN/k0_bcca)**(1.0/df_bcca)
    A0      = np.pi
    AC      = np.pi * rc * rc
    A_BCCA  = PN * A0 * minato06_bcca(PN)
    AC_BCCA = np.pi * rc_bcca * rc_bcca
    A       = 1.0/A_BCCA+1.0/AC-1.0/AC_BCCA
    A       = 1.0/A
    G       = min(A/PN/A0,1.0)

    return G