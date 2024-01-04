# mwprop.ne2001p v1.0 Jan 2024

from numpy import  log10

def calc_mwpsr_taudm(DM, A=2.98e-7, B=3.55e-5, alph=1.4, gam=3.1, sigma=0.76, sw=True):

    """
    Evaluates the fit to pulsar scattering

    tau(1 GHz, ms) = A \times DM^\alpha (1 + B \times DM^\beta) \times 10^{\pm sigma}

    Input: 
        DM = dispersion measure (pc/cc)
        A,B,alph,gam,sigma = values of fitted parameters
        sw = True  ==> use model as is = appropriate for Galactic objects
             False ==> multiply by xsw2pw = 3 = appropriate for extragalactic objects

    Output:
        taudm  = scattering time at 1 GHz in ms
        taudm_minus,plus = scattering time at \mp 1 sigma in log10
        taudm_outer_extend_minus,plus = extension of the low DM part

    """

    # xsw2pw = correction from spherical-wave scattering to plane-wave scattering
    xsw2pw = 3

    taudm = A * DM**alph * (1. + B * DM**gam)
    taudm_minus = 10.**(log10(taudm) - sigma)
    taudm_plus = 10.**(log10(taudm) + sigma)

    tau_outer_extend_minus= 10.**(log10(A*DM**alph) - sigma)
    tau_outer_extend_plus = 10.**(log10(A*DM**alph) + sigma)

    if sw is False:

        taudm *= xsw2pw
        taudm_minus *= xsw2pw
        taudm_plus *= xsw2pw

        tau_outer_extend_minus *= xsw2pw
        tau_outer_extend_plus *= xsw2pw
     
    return taudm, taudm_minus, taudm_plus, tau_outer_extend_minus, tau_outer_extend_plus

  
