import numpy as np

def minato06_bpca(PN):
    """ 
    --------------------------------------------------------------------------------
       Equation (A.1) in Minato et al. (2006, A&A, 452, 701) (BPCA)
    --------------------------------------------------------------------------------
     PN : number of monomers
     G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2.
          Thus, "G" is non-dimensional quantity.
    -------------------------------------------------------------------------------- 
    """

    if PN > 16:
        pca = 4.27*PN**(-0.315)*np.exp(-1.74/PN**0.243)
    else:
        pca = 12.5*PN**(-0.315)*np.exp(-2.53/PN**0.0920)

    return pca

def minato06_bcca(PN):
    """ 
    --------------------------------------------------------------------------------
       Equation (A.2) in Minato et al. (2006, A&A, 452, 701)  (BCCA)
    --------------------------------------------------------------------------------
     PN : number of monomers
     G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2.
          Thus, "G" is non-dimensional quantity.
    -------------------------------------------------------------------------------- 
    """

    if PN > 16:
        cca = 0.352 + 0.566*PN**(-0.138)
    else:
        cca = 12.5*PN**(-0.315)*np.exp(-2.53/PN**0.0920)

    return cca

def chain_exact(PN):
    """ 
    --------------------------------------------------------------------------------
        Exaxt solution of linear chain clusters
    --------------------------------------------------------------------------------
     PN : number of monomers
     G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2.
          Thus, "G" is non-dimensional quantity.
    -------------------------------------------------------------------------------- 
    """

    G = 1.0/PN + (1.0-1.0/PN)*8.0/np.pi/3.0

    return G

def chararea(PN,df,k0):
    """
    --------------------------------------------------------------------------------
        Characteristic cross section
    --------------------------------------------------------------------------------
     PN : number of monomers
     df : fractal dimension
     k0 : fractal prefactor
     G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2.
          Thus, "G" is non-dimensional quantity.
    --------------------------------------------------------------------------------
    """

    G = 5.0 * PN ** (-1.0+2.0/df)/3.0/k0**(2.0/df)
    G = min(G,1.0)
    return G

def okuzumi09(PN,df,k0):
    """
    --------------------------------------------------------------------------------
        Equation (47) in Okuzumi et al. (2009)
    --------------------------------------------------------------------------------
     PN : number of monomers
     df : fractal dimension
     k0 : fractal prefactor
     G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2.
          Thus, "G" is non-dimensional quantity.
    --------------------------------------------------------------------------------
    """

    k0_bcca = 1.04
    df_bcca = 1.9
    G_BCCA  = minato06_bcca(PN)
    GC      = 5.0 * PN ** (-1.0+2.0/df)/k0**(2.0/df)/3.0
    GC_BCCA = 5.0 * PN ** (-1.0+2.0/df_bcca)/k0_bcca**(2.0/df_bcca)/3.0
    G       = 1.0/G_BCCA + 1.0/GC - 1.0/GC_BCCA
    G       = min(1.0/G,1.0)
    return G

def meakin88(PN,G):
    """
    --------------------------------------------------------------------------------
        Equation (4) in Meakin and Donn (1988)
    --------------------------------------------------------------------------------
     PN : number of monomers
     G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2.
          Thus, "G" is non-dimensional quantity.
    --------------------------------------------------------------------------------
    """
    
    sig0 = 0.25 * pi
    A    = 0.2403
    B    = 0.5172
    beta = 0.8465
    sig  = (A+B*N**(beta-1.0))
    G    = sig/sig0

    return G

def paszun09(PN,df,k0):
    """
    !--------------------------------------------------------------------------------
    ! Equation (11) in Paszun & Dominik (2009, A&A, 507, 1023)
    ! by assuming that the outermost radius equals to the characteristic radius.
    !--------------------------------------------------------------------------------
    ! Input  : PN : Number of monomer particles
    !          k0 : Fractal prefactor
    !          df : Fractal dimension
    ! Return : G  : The geocross section of the aggregate normalized 
    !               by a factor N*pi*r0^2. 
    !
    ! Note:
    ! By assuming the outermost radius of the aggregate is the same as its characteristic
    ! radius, we obtain
    !
    ! R_A     / PN^{1.33+0.3/df}  \^{1/3.3}
    ! --- =  | ------------------- |,
    ! R_0     \ 1.21 * k0^{0.3/D} /
    !
    ! where R_A is the area-equivalent radius.
    !--------------------------------------------------------------------------------
    """

    pow1 = 1.33+0.3/df
    pow2 = -0.3/df
    ra   = (PN**pow1*k0**pow2/1.21)**(1.0/3.3)
    G    = ra * ra / PN 

    return G

