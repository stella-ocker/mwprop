# mwprop.ne2001p v1.0 Jan 2024

"""
density_ne2001.py

Returns the local electron density, fluctuation parameters, weights
and clump/void flags for the position (x, y, z)

Relacement for Fortran subroutine density_2001 


Comments from Fortran code density.NE2001.f
c final version of NE2001
c returns densities, F parameters and weights of the various components
c mods:
c 28 July 02:
c   put in 'if' statements in subroutine density_2001 so that
c   function calls are not done if particular weights
c   (wg1, wg2, etc.) are set to zero in gal01.inp
c   This is (a) cleaner and (b) much more efficient if
c   the clump or void component is flagged out.

Python version JMC 2020 Jan 28

Changes:
    2022 Jan 01: added functions density_2001_smallscale_comps and density_2001_smooth_comps
                 so that different LoS grids can be used to integrate different contributions
                 (to speed up calculations)
"""

from mwprop.ne2001p.config_ne2001p import *

from mwprop.ne2001p.density_components import * 
from mwprop.ne2001p.ne_lism import *
from mwprop.ne2001p.neclumpN_NE2001_fast import *
from mwprop.ne2001p.nevoidN_NE2001 import *
from mwprop.ne2001p.ne_arms_ne2001p import *

script_path = os.path.dirname(os.path.realpath(__file__))

# ---------------------------------------------------------------------------

def density_2001(x, y, z, inds_relevant=None, verbose=False):  
    """
    (Edited) Comments from Fortran code:

     Returns seven components of the free electron density of the
     interstellar medium at Galactic location (x,y,z).

     input:
      x, y, z = galactocentric location (kpc)
          Right-handed coordinate system
          x is in l=90 direction
          y is in l=180 direction
          The sun is at (x,y,z) = (0,rsun,0)
     output:
      electron densities in cm^{-3}:
      ne1:    outer, thick disk
      ne2:    inner, thin disk (annular in form)
      nea:    spiral arms
      negc:   galactic center component
          nelism: local ISM component
          necN:   contribution from discrete 'clumps'
          nevN:   contribution from voids
       fluctuation parameters (one for each ne component):
          F1, F2, Fa, Fgc, Flism, FcN, FvN
       flags:
          whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
             wlism: 1 if x,y,z is in any of the four LISM components
              wLDR: 1 if in LDR, 0 if not
              wLHB: 1 if in LHB, 0 if not
              wLSB: 1 if in LSB, 0 if not
            wLOOPI: 1 if in LoopI, 0 if not
          (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
          hitclump: clump number that x,y,z is in (0 if none)
           hitvoid: void number that x,y,z is in (0 if none)
    25 May 2002
    based on routines from TC93 and test routines from 1999-2002 by JMC.
    """

    # Assumes model parameters have been put into global dictionaries
    # Initiate values in case components are flagged out

    ne1 = ne2 = nea = negc = nelism = necN = nevN = 0
    wlism=wldr=wlhb=wlsb=wloopI = hitclump = hitvoid = wvoid = whicharm = 0


    if wg1 == 1: 
        ne1, F1 = ne_outer(x,y,z)
    if wg2 == 1: 
        ne2, F2 = ne_inner(x,y,z)
    if wga == 1: 
        nea, Fa, whicharm = ne_arms_ne2001p(x,y,z, Ncoarse=Ncoarse)
    else:
        nea = Fa = 0.
    if wggc == 1:  
        negc, Fgc = ne_gc(x,y,z) 
    if wglism == 1: 
        nelism, Flism, wlism, wldr, wlhb, wlsb, wloopI = ne_lism(x,y,z)
    if wgcN == 1: 
        necN, FcN, hitclump = neclumpN(x,y,z, inds_relevant=inds_relevant)
    if wgvN == 1: 
        nevN, FvN, hitvoid, wvoid = nevoidN(x,y,z)

    if verbose:
        print('density: ', ne1, ne2, F1, F2)

    
    return ne1,ne2,nea,negc,nelism,necN,nevN, \
           F1, F2, Fa, Fgc, Flism, FcN, FvN, \
           whicharm, wlism, wldr, wlhb, wlsb, wloopI, \
           hitclump, hitvoid, wvoid


# ---------------------------------------------------------------------------

def density_2001_smallscale_comps(x, y, z, inds_relevant=None, verbose=False):  
    """
    (Edited) Comments from Fortran code:

     Returns four components of the free electron density of the
     interstellar medium at Galactic location (x,y,z).

     input:
      x, y, z = galactocentric location (kpc)
          Right-handed coordinate system
          x is in l=90 direction
          y is in l=180 direction
          The sun is at (x,y,z) = (0,rsun,0)
     output:
      electron densities in cm^{-3}:
      negc:   galactic center component
      nelism: local ISM component
      necN:   contribution from discrete 'clumps'
      nevN:   contribution from voids
      fluctuation parameters (one for each ne component):
          Fgc, Flism, FcN, FvN
      flags:
             wlism: 1 if x,y,z is in any of the four LISM components
              wLDR: 1 if in LDR, 0 if not
              wLHB: 1 if in LHB, 0 if not
              wLSB: 1 if in LSB, 0 if not
            wLOOPI: 1 if in LoopI, 0 if not
          (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
          hitclump: clump number that x,y,z is in (0 if none)
           hitvoid: void number that x,y,z is in (0 if none)
    01 Jan 2022  based on fortran core from 25 May 2002
    """

    # Assumes model parameters have been put into global dictionaries
    # Initiate values in case components are flagged out

    negc = nelism = necN = nevN = 0
    wlism=wldr=wlhb=wlsb=wloopI = hitclump = hitvoid = wvoid = 0

    if wggc == 1:  
        negc, Fgc = ne_gc(x,y,z) 
    if wglism == 1: 
        nelism, Flism, wlism, wldr, wlhb, wlsb, wloopI = ne_lism(x,y,z)
    if wgcN == 1: 
        necN, FcN, hitclump, arg = neclumpN(x,y,z, inds_relevant=inds_relevant) # SKO added arg to output for debugging
    if wgvN == 1: 
        nevN, FvN, hitvoid, wvoid = nevoidN(x,y,z)

    if verbose:
        print('density: ', negc, nelism, necN, nevN) 

    if wgcN == 0:
        necN = 0
        FcN = 0
        hitclump = 0
    
    return negc,nelism,necN,nevN, \
           Fgc, Flism, FcN, FvN, \
           wlism, wldr, wlhb, wlsb, wloopI, \
           hitclump, hitvoid, wvoid

# ---------------------------------------------------------------------------

def density_2001_smooth_comps(x, y, z, verbose=False):  
    """
    Returns only the electron density for the smooth components.
    Utility: can be coarsely sampled to speed up computations.

    Output density is the sum ne1 + ne2 + nea at (x,y,z)

     input:
      x, y, z = galactocentric location (kpc)
          Right-handed coordinate system
          x is in l=90 direction
          y is in l=180 direction
          The sun is at (x,y,z) = (0,rsun,0)
     output:
      electron densities in cm^{-3}:
      ne1:    outer, thick disk
      ne2:    inner, thin disk (annular in form)
      nea:    spiral arms
    01 Jan 2022
    based on routines from NE2001 
    """

    # Assumes model parameters have been put into global dictionaries
    # Initiate values in case components are flagged out

    ne1 = ne2 = nea = 0
    F1 = F2 = Fa = 0
    whicharm = 0

    if wg1 == 1: 
        ne1, F1 = ne_outer(x,y,z)
    if wg2 == 1: 
        ne2, F2 = ne_inner(x,y,z)
    if wga == 1: 
        nea, Fa, whicharm = ne_arms_ne2001p(x,y,z, Ncoarse=Ncoarse)

    if verbose:
        print('density: ', ne1, ne2, F1, F2)

    # Define smooth components
    # Fsmooth defined so that 
    # Fsmooth*ne_smooth**2 = F1*ne1**2 + F2*ne2**2 + Fa*nea**2`
    wne1 = wg1*ne1
    wne2 = wg2*ne2
    wnea = wga*nea

    ne_smooth = wne1 + wne2 + wnea
    
    if ne_smooth > 0.:
        Fsmooth = (F1*(wne1)**2 + F2*wne2**2 + Fa*wnea**2) / ne_smooth**2
    else:
        Fsmooth = 0.

    return ne1,ne2,nea, F1, F2, Fa, whicharm, ne_smooth, Fsmooth

    """
    if verbose:
        return ne1,ne2,nea, \
               F1, F2, Fa, \
               whicharm, ne_smooth
    else:
        return ne_smooth
    """
