# mwprop.ne2001p v1.0 Jan 2024

"""
dmdsm_ne2001p.py
NE2001 integrator

Line of sight integration routines for NE2001p. 

Separate routines for calculating distance D from input dispersion measure DM
and vice verse.

This version is for NE2001p = Python version of NE2001 (fortran).
JMC 2021 Dec 28 - 2022 Jan 05
SKO 2022 Mar - 2023 Nov
JMC 2023 Dec 31 

This version includes:
    1. Use of coarse sampling + spline interpolation of smooth density components 
       (ne1, ne2, nea from thin and thick disks and spiral arm components)
    2. Fine sampling on all other components (LISM, GC, clumps, and voids)
    3. Ancillary analyses of the line of sight are done only if do_analysis = True
       in calling argument of dmdsm_dm2d in order to minimize computation times. .

Based on dmdsm_dm2d_ne2001_dev.py

Note mwprop needs to be in PYTHONPATH
"""
#import os
#script_path = os.path.dirname(os.path.realpath(__file__))

import numpy as np
from numpy import log10
from numpy import array, linspace, where, size
from numpy import digitize, interp, trapz

from mwprop.ne2001p.config_ne2001p import *
from mwprop.ne2001p.density_ne2001 import *
from mwprop.ne2001p.neclumpN_NE2001_fast import *

import datetime as datetime
import time

import sys,os

script_path = os.path.dirname(os.path.realpath(__file__))
basename = sys.argv[0].split('/')[-1].split('.')[0]
now = datetime.datetime.now()
plotstamp = basename + '_' + str(now).split('.')[0]

# ----------------------------------------------------------------------

def calc_galcentric_vecs(l, b, dmax, Ns):
    """
    Calculates vectors for line of sight and galactocentric coordinates
    Input:
        l,b = Galactic longitude and latitude (rad)
        dmax = maximum distance along LoS
        Ns = number of steps along LoS
    Output:
        svec = vector of steps along the LoS up to dmax
        xvec, yvec, zvec = vectors of x,y,z coordinates
    """
    sl=sin(l)
    cl=cos(l)
    sb=sin(b)
    cb=cos(b)

    svec = linspace(0, dmax, Ns)
    rvec = svec * cb
    xvec = rvec * sl
    yvec = rsun - rvec * cl
    zvec = svec * sb
    return svec, xvec, yvec, zvec

# ----------------------------------------------------------------------
# SKO -- 3/2/22 changed ds_coarse to 0.1 for both functions (was previously 0.2 for dm2d), changed ds_fine to 0.005
# changed Nsmin to 20 from 10
# SKO -- 3/6/22 -- changed ds_fine back to 0.01, issue was actually in nevoidN_NE2001.py

def dmdsm_dm2d(l, b, dm_target, ds_coarse=0.1, ds_fine=0.01, Nsmin=20, 
    dm2d_only = False, do_analysis=True, plotting=False, verbose=False, debug=False):
    """
    Integrates electron density from NE2001 model to reach the target DM.
    in the direction expressed in Galactic coordinates. 

    Computes pulsar distance and scattering measure
    from model of Galactic electron distribution.

    Input: l         galactic longitude in radians
           b         galactic latitude in radians
           dm_target input value DM (dispersion measure in pc/cm^3)
           ds_coarse coarse step size along LoS (kpc)
           ds_fine   fine step size along LoS (kpc)
           Nsmin     minimum number of samples to use along LoS
           dm2d_only True => calculate only the distance; otherwise also calculate SM.
           do_analysis True => analyze components of line of sight: dm, sm, lism, arms
           verbose   prints a few diagnostics

    Output:
           limit  (set to '>' if only a lower distance limit can be
               given; otherwise set to ' ')
           dist          calculated distance or input distance
           dmpsr         calculated DM or input DM
           sm            scattering measure, uniform weighting) (kpc/m^{20/3} 
           smtau         scattering measure, weighting for pulse broadening 
           smtheta       scattering measure, weighting for angular broadening
                          of galactic sources 
           smiso          scattering measure appropriate for calculating the
                           isoplanatic angle at the source's location
           Uses Kolmogorov spectral index = 11/3
           Useful constants: c_sm = (sikol - 3) / (2 * (2*pi)**(4-sikol)) 
                             sm_factor = c_sm * units_conversion to kpc m^{-20/3} 
                             (defined in config_ne2001p.py)       
    """
    sm_iso_index = sikol - 2                # should be 5/3 for Kolmogorov index sikol=11/3
    limit=' '                       

    # Do multiple passes on integration:
    # Pass  0: coarse integration to get dmax to use for fine steps
    # Pass >0: distance estimate using fine distance steps

    dm_reached = False

    npasses = 3
    mult_dmaxnom = 1.5          # multiplier of dhat to guard against undershooting

    npass = 0
    dhat = 0
    dm_calc_max = 0
    while  npass < npasses and dm_reached is False:
 
        if npass==0:
            # Calculate nominal maximum distance to calculate initial set of samples
            # Use small nominal n_e for small DM, larger for large DM
            # Rationale:  large DM requires inner Galaxy LoS that sample larger n_e
            ne_nom = 0.01 * (1 + log10(dm_target))
            dmax_integrate = \
                 min(mult_dmaxnom * dm_target/ ne_nom / pc_in_kpc, dmax_ne2001p_integrate)
            Ns_coarse = int(dmax_integrate/ds_coarse)
            if Ns_coarse < Nsmin:
                Ns_coarse = Nsmin
                ds_coarse = dmax_integrate/ Ns_coarse
            Ns_fine = int(dmax_integrate/ds_fine)
            if Ns_fine < Nsmin:
                Ns_fine = Nsmin
                ds_fine = dmax_integrate / Ns_fine
                
            # Identify clumps that contribute to this LoS:
            # SKO 11/23 -- rcmult now defined in config_ne2001.py; rcmult smaller for clumps with edge = 0 (hard cutoff) 
            relevant_clump_indices = relevant_clumps(l, b, dmax_integrate, rcmult)
           
            if debug:
                print(np.size(relevant_clump_indices), 
                    ' relevant clumps out of ', nclumps, ' total')
                print(relevant_clump_indices)
        else:
            dmax_integrate = min(mult_dmaxnom * dhat, dmax_ne2001p_integrate)
            Ns_fine = int(dmax_integrate / ds_fine) 
            if Ns_fine < Nsmin:
                Ns_fine = Nsmin 
                ds_fine = dmax_integrate / Ns_fine
        if debug:
            print('dmax_integrate = ', dmax_integrate)
            print('Ns_fine, ds_fine = ', Ns_fine, ds_fine)
            print('ds_coarse = ', ds_coarse, ' ds_fine = ', ds_fine)
            print('Ns_coarse = ', Ns_coarse, ' Ns_fine = ', Ns_fine)
 
        sc_vec, xc_vec, yc_vec, zc_vec = calc_galcentric_vecs(l, b, dmax_integrate, Ns_coarse)
        sf_vec, xf_vec, yf_vec, zf_vec = calc_galcentric_vecs(l, b, dmax_integrate, Ns_fine)

        # Obtain electron density components, F parameters, weights, etc.
        # Use separate calls to density_2001_smooth_comps and density_2001_smallscale_comps

        # ----------------------------------------------
        # Smooth, large-scale components on coarse grid: 
        # ----------------------------------------------
        # Note only cnd_smooth, cFsmooth needed here
        cne1,cne2,cnea, cF1, cF2, cFa, cwhicharm, cne_smooth, cFsmooth = \
            array([
               density_2001_smooth_comps(xc_vec[j],yc_vec[j],zc_vec[j]) 
               for j in range(Ns_coarse)
            ]).T

        # Spline functions:
        cs_ne_smooth = CubicSpline(sc_vec, cne_smooth)
        cs_F_smooth = CubicSpline(sc_vec, cFsmooth)

        # Evaluate smooth components on fine grid  using spline functions:
        ne_smooth = cs_ne_smooth(sf_vec)
        F_smooth = cs_F_smooth(sf_vec)

        # Resample cwhicharm using digitize: use coarse vec as bins for fine vec:
        inds_whicharm = np.digitize(sf_vec, sc_vec, right=True)    
        whicharm = cwhicharm[inds_whicharm] 

        # ------------------------------------
        # Small-scale components on fine grid:
        # ------------------------------------
        negc,nelism,necN,nevN, Fgc, Flism, FcN, FvN, wlism, wldr, wlhb, wlsb, wloopI, \
           hitclump, hitvoid, wvoid = \
               array([
                   density_2001_smallscale_comps(\
                       xf_vec[j],yf_vec[j],zf_vec[j], inds_relevant=relevant_clump_indices) \
                       for j in range(Ns_fine)\
               ]).T

        wtotal = (1-wgvN*wvoid)*(1-wglism*wlism)        # used for SM calculations
        ne_ex_clumps_voids = (1.-wglism*wlism) * (ne_smooth + wggc*negc) + wglism*wlism*nelism
        ne = (1-wgvN*wvoid)*ne_ex_clumps_voids  + wgvN*wvoid*nevN + wgcN*necN 

        dm_cumulate_vec = \
            pc_in_kpc * array([trapz(ne[:j], sf_vec[:j]) for j in range(1, Ns_fine+1) ])
        dm_calc_max = dm_cumulate_vec[-1]       # maximum dm calculated for this pass

        # Interpolate to get distance estimate:
        dhat = interp(dm_target, dm_cumulate_vec, sf_vec)

        if debug: 
            print('Pass ', npass, '  Ns_fine = %d  dmax_integrate = %5.2f   dhat = %5.2f'\
                %(Ns_fine, dmax_integrate, dhat))
        if dm_calc_max < dm_target: # SKO -- modified warning
            #print('Warning: Pass %d    terminal DM < target DM (%6.1f, %6.1f)   dhat = %6.2f' %(npass, dm_calc_max, dm_target, dhat)) 
            if npass+1 == npasses:
                print('Warning: terminal DM < target DM')
            dhat *= mult_dmaxnom                # increase dhat to reach target DM
        else:
            dm_reached = True
        npass += 1

    # Convert flags to integers:
    whicharm = whicharm.astype(int)
    hitclump = hitclump.astype(int)
    hitvoid = hitvoid.astype(int)

    # Integrate using trapz to get cumulative DM:
    dm_cumulate_vec = pc_in_kpc * array([trapz(ne[:j], sf_vec[:j]) for j in range(1,Ns_fine+1)])
    dm_calc_max = dm_cumulate_vec[-1]       # maximum dm calculated 

    # Test if calculated DMs reach dm_target.
    # If not, attribute this to a lower bound on the distance, as in NE2001 (fortran),
    # though it might be better to attribute this to insufficent electrons.
    # In new version NE2001x, use this second approach. 

    if dm_calc_max < dm_target:     # did not reach target DM; find dhat lower limit
     
        limit = '>'
        indlim = where(dm_cumulate_vec==dm_cumulate_vec[-1])[0][0]
        dm_reached = dm_cumulate_vec[indlim]
        dhat_lower_limit = sf_vec[indlim]
        dhat = dhat_lower_limit     # for consistency with return statements
        dm_return = dm_reached
        dhat_return = dhat_lower_limit
       
        #if debug:
            #indmin = where(dm_cumulate_vec==dm_cumulate_vec[-1])[0][0]
            #print('indmin = ', indmin, sf_vec[indmin])
    else:                           # Reached target DM; interpolate to get distance estimate:
        dhat = interp(dm_target, dm_cumulate_vec, sf_vec)
        dhat_return = dhat
        dm_return = dm_target


    if plotting:
        plot_dm_along_LoS(
            dm_target, dhat, sf_vec, dm_cumulate_vec, relevant_clump_indices, dc, 
            plot_dm_target=True, saveplot=True)

    if dm2d_only is True:
        if do_analysis:

            savedir = os.getcwd()+'/output_ne2001p/'
            f24outfile  = savedir+'f24_dm2d_only_lism_etc_paths' + '.txt'
            f24 = open(f24outfile, 'w')

            f25outfile  = savedir+'f25_dm2d_only_ne_vs_s' + '.txt'
            f25 = open(f25outfile, 'w')

            sf_vec, ne, ne1, ne2, nea = \
                analysis_dmd_dm_only(f24, f25,
                    l, b, dm_target, dhat, sf_vec, xf_vec, yf_vec, zf_vec, sc_vec,
                    cne1, cne2, cnea, ne, whicharm, hitclump, hitvoid, dm_cumulate_vec,
                    nelism, negc, necN, nevN, wtotal, wlism, wvoid,
                    wlhb, wldr, wlsb, wloopI)
        if debug:
            return limit, dhat_return, dm_return, sf_vec, dm_cumulate_vec
        else:
            return limit, dhat_return, dm_return

    else:           # then also calculate scattering measure, do_analysis, etc.

        # Calculate dSM (proportional to Cn2) from individual dSM components
        # Note: dsm1, dsm2, dsma are subsumed by dsm_smooth
        # dsm1, dsm2, dsma are calculated in the `analysis' code block (if requested)

        # For now, the dsm quantities are not multiplied by sm_factor.
        # That is done later in calculating sm quantities from dsm quantities.
        # Could change this to make it a bit more transparent
        # e.g. multiply dsm quantities by sm_factor here.

        dsm_smooth = wtotal * ne_smooth**2 * F_smooth
        dsmgc = wtotal*wggc*negc**2*Fgc
        dsmlism = (1.-wgvN*wvoid)*wglism*wlism*nelism**2*Flism
        dsmcN = wgcN*necN**2*FcN
        dsmvN = wgvN*wvoid*nevN**2*FvN
        dsm = dsm_smooth + dsmgc + dsmlism + dsmcN + dsmvN

        # Calculate integrals needed to evaluate SM.
        # Integrate from s = 0 to s = dhat (starting from observer's position at s = 0). 
        # First integrate over sf_vec then use cubic spline to find SM at d = dhat.

        Nsf1 = np.size(sf_vec) + 1
        dsm_cumulate1_vec = array([trapz(dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])
        dsm_cumulate2_vec = \
            array([trapz(sf_vec[:j]*dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])
        dsm_cumulate3_vec = \
            array([trapz(sf_vec[:j]**2*dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])
        dsm_cumulate4_vec = \
            array([trapz(sf_vec[:j]**sm_iso_index*dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])

        sm_cumulate = sm_factor * dsm_cumulate1_vec
        smtau_cumulate = 6 * (sm_factor/dhat) * (dsm_cumulate2_vec - dsm_cumulate3_vec/dhat)
        smtheta_cumulate = 3 * sm_factor * \
            (dsm_cumulate1_vec + dsm_cumulate3_vec/dhat**2 - 2*dsm_cumulate2_vec/dhat)
        smiso_cumulate = sm_factor * dsm_cumulate4_vec

        # cubic splines:

        sm_cs = CubicSpline(sf_vec, sm_cumulate)
        smtau_cs = CubicSpline(sf_vec, smtau_cumulate)
        smtheta_cs = CubicSpline(sf_vec, smtheta_cumulate)
        smiso_cs = CubicSpline(sf_vec, smiso_cumulate)

        sm_hat = sm_cs(dhat)
        smtau_hat = smtau_cs(dhat)
        smtheta_hat = smtheta_cs(dhat)
        smiso_hat = smiso_cs(dhat)

        if plotting:
            plot_dm_sm_along_LoS(
                dm_target, dm_cumulate_vec, dhat, sf_vec, 
                sm_hat, smtau_hat, smtheta_hat, smiso_hat,
                dsm, sm_cumulate, smtau_cumulate, smtheta_cumulate, smiso_cumulate,
                saveplot=True)

        if do_analysis:

            savedir = os.getcwd()+'/output_ne2001p/'
            f24outfile  = savedir+'f24_dm2d_lism_etc_paths' + '.txt'
            f24 = open(f24outfile, 'w')

            f25outfile  = savedir+'f25_dm2d_ne_dsm_vs_s' + '.txt'
            f25 = open(f25outfile, 'w')

            svec, ne, ne1, ne2, nea = \
            analysis_dmd_dm_and_sm(f24, f25,
                l, b, dm_target, dhat, sf_vec, xf_vec, yf_vec, zf_vec, sc_vec,
                cne1, cF1, cne2, cF2, cnea, cFa,
                sm_hat, ne, dsm, whicharm, hitclump, hitvoid, dm_cumulate_vec,
                nelism, negc, necN, nevN, wtotal, wlism, wvoid,
                dsmgc, dsmlism, dsmcN, dsmvN,
                wlhb, wldr, wlsb, wloopI)

        if debug:
            return limit, dhat_return, dm_return, sm_hat, smtau_hat, smtheta_hat, smiso_hat, \
               sf_vec, dm_cumulate_vec
        else:
            return limit, dhat_return, dm_return, sm_hat, smtau_hat, smtheta_hat, smiso_hat

# ----------------------------------------------------------------------

def dmdsm_d2dm(l, b, d_target, ds_coarse, ds_fine, Nsmin, 
    d2dm_only=False, do_analysis=True, plotting=False, verbose=False): # ds_coarse, ds_fine, Nsmin defined in config_ne2001

    """
    Integrates electron density from NE2001 model out to a specified
    distance  in the direction expressed in Galactic coordinates. 

    Computes pulsar distance and scattering measure
    from model of Galactic electron distribution.

    Input: l         galactic longitude in radians
           b         galactic latitude in radians
           d_target  input distance to integrate to (kpc)
           ds_coarse coarse step size along LoS (kpc)
           ds_fine   fine step size along LoS (kpc)
           Nsmin     minimum number of samples to use along LoS
           d2dm_only True => calculate DM only; otherwise also calculate SM.
           do_analysis True => analyze components of line of sight: dm, lism, arms (sm)
           plotting  True => plot dm etc vs distance along path
           verbose   N/A presently

    Output:
           dist          calculated distance or input distance
           dm            calculated DM 
           sm            scattering measure, uniform weighting) (kpc/m^{20/3} 
           smtau         scattering measure, weighting for pulse broadening 
           smtheta       scattering measure, weighting for angular broadening
                          of galactic sources 
           smiso          scattering measure appropriate for calculating the
                           isoplanatic angle at the source's location
           Uses Kolmogorov spectral index = 11/3
           Useful constants: c_sm = (sikol - 3) / (2 * (2*pi)**(4-sikol)) 
                             sm_factor = c_sm * units_conversion to kpc m^{-20/3} 
                             (defined in config_ne2001p.py)       
    """
    # ----------------------------------------------------------------------

    sm_iso_index = sikol - 2                # should be 5/3 for Kolmogorov index sikol=11/3
    limit=' '                               # placeholder for legacy code

    # Coarse and fine integration grids 
    Ns_coarse = int(d_target / ds_coarse)
    if Ns_coarse < Nsmin:
        Ns_coarse = Nsmin
        ds_coarse = d_target / Ns_coarse

    Ns_fine = int(d_target / ds_fine)
    if Ns_fine < Nsmin:
        Ns_fine = Nsmin
        ds_fine = d_target / Ns_fine

    relevant_clump_indices = relevant_clumps(l, b, d_target, rcmult)

    sc_vec, xc_vec, yc_vec, zc_vec = calc_galcentric_vecs(l, b, d_target, Ns_coarse)
    sf_vec, xf_vec, yf_vec, zf_vec = calc_galcentric_vecs(l, b, d_target, Ns_fine)

    # ---------------------------------------------
    # Smooth, large-scale components on coarse grid
    # ---------------------------------------------
    # Note only cnd_smooth, cFsmooth needed here:
    cne1,cne2,cnea, cF1, cF2, cFa, cwhicharm, cne_smooth, cFsmooth = \
        array([
           density_2001_smooth_comps(xc_vec[j],yc_vec[j],zc_vec[j]) for j in range(Ns_coarse)
        ]).T

    # Spline functions:
    cs_ne_smooth = CubicSpline(sc_vec, cne_smooth)
    cs_F_smooth = CubicSpline(sc_vec, cFsmooth)

    # Evaluate using spline functions:
    ne_smooth = cs_ne_smooth(sf_vec)
    F_smooth = cs_F_smooth(sf_vec)

    # Resample cwhicharm using digitize: use coarse vec as bins for fine vec:
    inds_whicharm = np.digitize(sf_vec, sc_vec, right=True)
    whicharm = cwhicharm[inds_whicharm]

    # ------------------------------------
    # Small-scale components on fine grid:
    # ------------------------------------
    # SKO -- tested changing inds_relevant to None from relevant_clump_indices
    negc,nelism,necN,nevN, Fgc, Flism, FcN, FvN, wlism, wldr, wlhb, wlsb, wloopI, \
       hitclump, hitvoid, wvoid = \
           array([
               density_2001_smallscale_comps(\
                   xf_vec[j],yf_vec[j],zf_vec[j], inds_relevant=relevant_clump_indices) \
                   for j in range(Ns_fine)\
           ]).T
    
    wtotal = (1-wgvN*wvoid)*(1-wglism*wlism)        # used for SM calculations
    ne_ex_clumps_voids = (1.-wglism*wlism) * (ne_smooth + wggc*negc) + wglism*wlism*nelism
    ne = (1-wgvN*wvoid)*ne_ex_clumps_voids  + wgvN*wvoid*nevN + wgcN*necN

    dm_cumulate_vec = \
        pc_in_kpc * array([trapz(ne[:j], sf_vec[:j]) for j in range(1, Ns_fine+1) ])
    dm_calc_max = dm_cumulate_vec[-1]       

    # floats -> ints:
    whicharm = whicharm.astype(int)
    hitclump = hitclump.astype(int)
    hitvoid = hitvoid.astype(int)

    if plotting:
        plot_dm_along_LoS(
            dm_calc_max, d_target, sf_vec, dm_cumulate_vec, relevant_clump_indices, dc, 
            which='d2dm', plot_dm_target=True, saveplot=True)

    if d2dm_only is True:
            if do_analysis:
                savedir = os.getcwd()+'/output_ne2001p/'
                f24outfile  = savedir+'f24_d2dm_only_lism_etc_paths' + '.txt'
                f24 = open(f24outfile, 'w')

                f25outfile  = savedir+'f25_d2dm_only_ne_vs_s' + '.txt'
                f25 = open(f25outfile, 'w')

                dhat = d_target             # so that code here is same as in dmdsm_dm2d
                dm_target = dm_calc_max     # ditto
                svec, ne, ne1, ne2, nea = \
                    analysis_dmd_dm_only(f24, f25,
                    l, b, dm_target, dhat, sf_vec, xf_vec, yf_vec, zf_vec, sc_vec,
                    cne1, cne2, cnea, ne, whicharm, hitclump, hitvoid, dm_cumulate_vec,
                    nelism, negc, necN, nevN, wtotal, wlism, wvoid,
                    wlhb, wldr, wlsb, wloopI)
            return limit, d_target, dm_calc_max
    else:

        dhat = d_target             # so that code here is same as in dmdsm_dm2d
        dm_target = dm_calc_max     # ditto

        # Calculate dSM (proportional to Cn2) from individual dSM components
        # Note: dsm1, dsm2, dsma are subsumed by dsm_smooth
        # dsm1, dsm2, dsma are calculated in the `analysis' code block (if requested)

        # For now, the dsm quantities are not multiplied by sm_factor.
        # That is done later in calculating sm quantities from dsm quantities.
        # Could change this to make it a bit more transparent 
        # e.g. multiply dsm quantities by sm_factor here.

        dsm_smooth = wtotal * ne_smooth**2 * F_smooth
        dsmgc = wtotal*wggc*negc**2*Fgc
        dsmlism = (1.-wgvN*wvoid)*wglism*wlism*nelism**2*Flism
        dsmcN = wgcN*necN**2*FcN
        dsmvN = wgvN*wvoid*nevN**2*FvN
        dsm = dsm_smooth + dsmgc + dsmlism + dsmcN + dsmvN

        # Calculate integrals needed to evaluate SM;
        # Note integrals start at observer's position
        # Integrate from s = 0 to s = dhat
        # Accomplish this by integrating over full svec used
        # and then use cubic spline to find SM at d = dhat:.

        Nsf1 = np.size(sf_vec) + 1
        dsm_cumulate1_vec = array([trapz(dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])
        dsm_cumulate2_vec = \
            array([trapz(sf_vec[:j]*dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])
        dsm_cumulate3_vec = \
            array([trapz(sf_vec[:j]**2*dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])
        dsm_cumulate4_vec = \
            array([trapz(sf_vec[:j]**sm_iso_index*dsm[:j], sf_vec[:j]) for j in range(1, Nsf1)])

        # sm quantities have proper SM units:

        sm_cumulate = sm_factor *dsm_cumulate1_vec
        smtau_cumulate = 6 * (sm_factor/dhat) * (dsm_cumulate2_vec - dsm_cumulate3_vec/dhat)
        smtheta_cumulate = 3 * sm_factor * \
            (dsm_cumulate1_vec + dsm_cumulate3_vec/dhat**2 - 2*dsm_cumulate2_vec/dhat)
        smiso_cumulate = sm_factor * dsm_cumulate4_vec

        # cubic splines:

        sm_cs = CubicSpline(sf_vec, sm_cumulate)
        smtau_cs = CubicSpline(sf_vec, smtau_cumulate)
        smtheta_cs = CubicSpline(sf_vec, smtheta_cumulate)
        smiso_cs = CubicSpline(sf_vec, smiso_cumulate)

        sm_hat = sm_cs(dhat)
        smtau_hat = smtau_cs(dhat)
        smtheta_hat = smtheta_cs(dhat)
        smiso_hat = smiso_cs(dhat)

        if plotting:
            plot_dm_sm_along_LoS(
                dm_target, dm_cumulate_vec, dhat, sf_vec,
                sm_hat, smtau_hat, smtheta_hat, smiso_hat,
                dsm, sm_cumulate, smtau_cumulate, smtheta_cumulate, smiso_cumulate,
                which='d2dm', saveplot=True)

        if do_analysis:

            savedir = os.getcwd()+'/output_ne2001p/'

            f24outfile  = savedir+'f24_d2dm_lism_etc_paths' + '.txt'
            f24 = open(f24outfile, 'w')

            f25outfile  = savedir+'f25_d2dm_ne_dsm_vs_s' + '.txt'
            f25 = open(f25outfile, 'w')

            svec, ne, ne1, ne2, nea = \
            analysis_dmd_dm_and_sm(f24, f25, 
                l, b, dm_target, dhat, sf_vec, xf_vec, yf_vec, zf_vec, sc_vec,
                cne1, cF1, cne2, cF2, cnea, cFa,
                sm_hat, ne, dsm, whicharm, hitclump, hitvoid, dm_cumulate_vec,
                nelism, negc, necN, nevN, wtotal, wlism, wvoid,
                dsmgc, dsmlism, dsmcN, dsmvN,
                wlhb, wldr, wlsb, wloopI)

        return limit, dhat, dm_target, sm_hat, smtau_hat, smtheta_hat, smiso_hat

# ----------------------------------------------------------------------

def analysis_dmd_dm_only(f24, f25, 
    l, b, dm_target, dhat, sf_vec, xf_vec, yf_vec, zf_vec, sc_vec,
    cne1, cne2, cnea, ne, whicharm, hitclump, hitvoid, dm_cumulate_vec,
    nelism, negc, necN, nevN, wtotal, wlism, wvoid,
    wlhb, wldr, wlsb, wloopI):
    
        """
        Assess contributions from different model components
        and print to files. 
 
        This version excludes any scattering variables
        (a selection made for faster execution)

        Note that some quantities are defined only here because otherwise
        they would slow down the dm->d calculation.

        Output files are defined outside this function (f24, f25)
        Units:
            l,b                                            rad
            dm_target, dm_cumulate_vec                     pc/cc           
            dhat, sfvec, xfvec, yf_vec, zf_vec, sc_vec     kpc
            cne1, cne2, cnea, ne, nelism, negc, necN, nevN 1/cc
            cF1, cF2, cFa                                  pc^{-2/3}  check

        Dimensionless:
            whicharm, hitclump, hitvoid, wtotal, wlism, wvoid, wlhb, wldr, wlsb, wloopI

        """
        # Calculate smooth components on fine grid using CubicSpline to produce
        # spline functions from coarse samples:

        cs_ne1 = CubicSpline(sc_vec, cne1)
        cs_ne2 = CubicSpline(sc_vec, cne2)
        cs_nea = CubicSpline(sc_vec, cnea)

        # Evaluate using spline functions:
        ne1 = cs_ne1(sf_vec)
        ne2 = cs_ne2(sf_vec)
        nea = cs_nea(sf_vec)

        # print differential values along LoS to file:
        print('DM target: ', dm_target, file=f25)
        print(rad2deg(l), rad2deg(b), ' ldeg, bdeg', file=f25)

        # Entries with '---' are to have placeholders so that output number of fields
        # is the same as for analysis that includes SM.
        print('    d       x       y       z       ne     -------  w    c   v  t     dm       nea     ---', file=f25)

        testdm = where(dm_cumulate_vec <= dm_target, 0, 1)     # 0 if dm < dm_target

        Ns_fine = sf_vec.size
        # print out values only up to one step beyond where sf_vec = dhat:
        inds = where(testdm==0)
        Ns_print = np.min((np.size(sf_vec), np.size(inds)+1))
 
        for n in range(Ns_print):
                    print(
                '{:7.3f}{:8.3f}{:8.3f}{:8.3f}{:9.4f}{:11.4e}{:2d}{:5d}{:4d}{:3b}{:9.3f}{:10.4f}{:8.2f}'.format(sf_vec[n],xf_vec[n],yf_vec[n],zf_vec[n],ne[n],0.,whicharm[n],hitclump[n],hitvoid[n], testdm[n], dm_cumulate_vec[n], nea[n], 0.), file=f25)
        f25.close()

        # Print LoS summary analysis to file 24:

        # --------------------------------------
        # Calculate individual end-point values:
        # Uses CubicSplines
        # --------------------------------------

        ddm1 = wtotal*wg1*ne1 * 1000
        ddm2 = wtotal*wg2*ne2 * 1000
        ddma = wtotal*wga*nea * 1000
        ddmgc = wtotal*wggc*negc * 1000
        ddmlism = (1.-wgvN*wvoid)*wglism*wlism*nelism * 1000
        ddmcN = wgcN*necN * 1000
        ddmvN = wgvN*wvoid*nevN * 1000

        dm1run = array([trapz(ddm1[0:j], sf_vec[0:j]) for j in range(size(ddm1))])
        dm2run = array([trapz(ddm2[0:j], sf_vec[0:j]) for j in range(size(ddm2))])
        dmarun = array([trapz(ddma[0:j], sf_vec[0:j]) for j in range(size(ddma))])
        dmgcrun = array([trapz(ddmgc[0:j], sf_vec[0:j]) for j in range(size(ddmgc))])
        dmlismrun = array([trapz(ddmlism[0:j], sf_vec[0:j]) for j in range(size(ddmlism))])
        dmcNrun = array([trapz(ddmcN[0:j], sf_vec[0:j]) for j in range(size(ddmcN))])
        dmvNrun = array([trapz(ddmvN[0:j], sf_vec[0:j]) for j in range(size(ddmvN))])

        cs_dm1 = CubicSpline(sf_vec, dm1run)
        cs_dm2 = CubicSpline(sf_vec, dm2run)
        cs_dma = CubicSpline(sf_vec, dmarun)
        cs_dmgc = CubicSpline(sf_vec, dmgcrun)
        cs_dmlism = CubicSpline(sf_vec, dmlismrun)
        cs_dmcN = CubicSpline(sf_vec, dmcNrun)
        cs_dmvN = CubicSpline(sf_vec, dmvNrun)

        dm1 = cs_dm1(dhat)
        dm2 = cs_dm2(dhat)
        dma = cs_dma(dhat)
        dmgc = cs_dmgc(dhat)
        dmlism = cs_dmlism(dhat)
        dmcN = cs_dmcN(dhat)
        dmvN = cs_dmvN(dhat)


        #---------------
        # LISM analysis:
        #---------------
        # Calculate path lengths through LISM components:
        # Take into account the weighting hierarchy, LHB:LOOPI:LSB:LDR

        ldr_path = lhb_path = lsb_path = loopI_path = 0.
        ldr_dist = lhb_dist = lsb_dist = loopI_dist = 0.

        lhb_path = trapz(wlhb, sf_vec)
        loopI_path = trapz((1-wlhb) * wloopI, sf_vec)
        lsb_path = trapz((1-wlhb)*(1-wloopI) * wlsb, sf_vec)
        ldr_path = trapz((1-wlhb)*(1-wloopI)*(1-wlsb) * wldr, sf_vec)

        if lhb_path > 0:
            lhb_dist = trapz(wlhb*sf_vec, sf_vec) / lhb_path 
        if loopI_path > 0:
            loopI_dist = trapz((1-wlhb)*wloopI*sf_vec, sf_vec) / loopI_path 
        if lsb_path > 0:
            lsb_dist = trapz((1-wlhb)*(1-wloopI)*wlsb*sf_vec, sf_vec) / lsb_path
        if ldr_path > 0:
            ldr_dist = trapz((1-wlhb)*(1-wloopI)*(1-wlsb)*wldr*sf_vec, sf_vec) / ldr_path

        # --------------------
        # Spiral arm analysis:
        # --------------------
        # Calculate path lengths through spiral arms:
        """
        pathlengths: 
            whicharm = 0,5 (currently).
                       1,4 for the equivalent of the TC93 arms
                         5 for the local arm
                         0 means interarm path segments
        note narmsmax = number of arms
        armpaths has one additional entry to include interarm path length, so
        narmsmax1 = narmsmax + 1 (set in config_ne2001p.py)
        """
        armpaths = zeros(narmsmax1)
        armdistances = zeros(narmsmax1)

        # Using for loop for greater code clarity 
        # (could use list comprehensions but expressions are clunky)
        for j in range(narmsmax1):
            inds = where(whicharm==j) 
            sinds = sf_vec[inds]
            winds = whicharm[inds] / max(j, 1)              # avoid div by zero for j=0
            armpaths[j] = trapz(winds, sinds)
            armdistances[j] = trapz(sinds*winds, sinds) / where(armpaths[j] > 0, armpaths[j], 1) 

        # Printing to f24 file:
        print('LISM path lengths (kpc) with weighting hierarchy LHB:LOOPI:LSB:LDR',
        file = f24)
        print('{}'.format(18*' '+ '  LHB     LoopI     LSB      LDR'.expandtabs(15)),
        file = f24)
        print(3*' '+  '{0:15s}{1:6.3f}{2:9.3f}{3:9.3f}{4:9.3f}'.format('Length', lhb_path, loopI_path, lsb_path, ldr_path), file=f24)
        print(3*' '+  '{0:15s}{1:6.3f}{2:9.3f}{3:9.3f}{4:9.3f}'.format('Mean Dist.', lhb_dist, loopI_dist, lsb_dist, ldr_dist), file=f24)

        print('\nFractional contributions to DM:', file=f24)
        print('  outer   inner    arms     gc    lism    clumps  voids       DM', file= f24)
        print('{:7.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:11.3f}'.format(
            dm1/dm_target, dm2/dm_target, dma/dm_target, dmgc/dm_target,
            dmlism/dm_target, dmcN/dm_target, dmvN/dm_target, dm_target), file=f24)

        print('\nPath lengths through spiral arms:', file=f24)
        print('  Arm      Mean Distance       Path Length    (arm=0 => interarm)', file=f24)

        for j in range(narmsmax1):
           print('  {:2d}{:18.3f}{:18.3f}'.format(j, armdistances[j], armpaths[j]), file=f24)

        f24.close()
        return sf_vec, ne, ne1, ne2, nea


# ----------------------------------------------------------------------

def analysis_dmd_dm_and_sm(f24, f25,
    l, b, dm_target, dhat, sf_vec, xf_vec, yf_vec, zf_vec, sc_vec,
    cne1, cF1, cne2, cF2, cnea, cFa,
    sm, ne, dsm, whicharm, hitclump, hitvoid, dm_cumulate_vec,
    nelism, negc, necN, nevN, wtotal, wlism, wvoid,
    dsmgc, dsmlism, dsmcN, dsmvN,
    wlhb, wldr, wlsb, wloopI):
    
    """
    Assess contributions from different model components
    and print to files. 

    Note that some quantities are defined only here because otherwise
    they would slow down the dm->d calculation.

    e.g. dsm1, dsm2, dsma

    Output files are defined outside this function (f24, f25)

    Units:
        l,b                                            rad
        dm_target, dm_cumulate_vec                     pc/cc
        dhat, sfvec, xfvec, yf_vec, zf_vec, sc_vec     kpc
        cne1, cne2, cnea, ne, nelism, negc, necN, nevN 1/cc
        cF1, cF2, cFa                                  pc^{-2/3}  check
        sm                                             kpc m^{-20/3}
        dsm, dsmgc, dsmlism, dsmcN, dsmvN              (sm units) / sm_factor

        Dimensionless: 
            whicharm, hitclump, hitvoid, wtotal, wlism, wvoid, wlhb, wldr, wlsb, wloopI
    """

    # Calculate smooth components on fine grid using CubicSpline to produce
    # spline functions from coarse samples:

    cs_ne1 = CubicSpline(sc_vec, cne1)
    cs_F1 = CubicSpline(sc_vec, cF1)
    cs_ne2 = CubicSpline(sc_vec, cne2)
    cs_F2 = CubicSpline(sc_vec, cF2)
    cs_nea = CubicSpline(sc_vec, cnea)
    cs_Fa = CubicSpline(sc_vec, cFa)

    # Evaluate using spline functions:
    ne1 = cs_ne1(sf_vec)
    ne2 = cs_ne2(sf_vec)
    nea = cs_nea(sf_vec)
    F1 = cs_F1(sf_vec)
    F2 = cs_F2(sf_vec)
    Fa = cs_Fa(sf_vec)

    # print differential values along LoS to file:
    print('DM target: ', dm_target, file=f25)
    print(rad2deg(l), rad2deg(b), ' ldeg, bdeg', file=f25)

    print('    d       x       y       z       ne         Cn2     w    c   v  t     dm       nea      Fa', file=f25)

    testdm = where(dm_cumulate_vec <= dm_target, 0, 1)     # 0 if dm < dm_target

    Ns_fine = sf_vec.size
    # print out values only up to one beyond where sf_vec = dhat:
    inds = where(testdm==0)
    Ns_print = np.min((np.size(sf_vec), np.size(inds)+1))

    for n in range(Ns_print):
        print(
                '{:7.3f}{:8.3f}{:8.3f}{:8.3f}{:9.4f}{:14.4e}{:2d}{:5d}{:4d}{:3b}{:9.3f}{:10.4f}{:8.2f}'.format(sf_vec[n], xf_vec[n], yf_vec[n], zf_vec[n], ne[n], sm_factor*dsm[n], whicharm[n], hitclump[n], hitvoid[n], testdm[n], dm_cumulate_vec[n], nea[n], Fa[n]), file=f25)
    f25.close()

    # print LoS summary analysis to file

    # --------------------------------------
    # Calculate individual end-point values:
    # Uses CubicSplines
    # --------------------------------------

        # ddm/ds components at fine resolution:
    ddm1 = wtotal*wg1*ne1 * 1000
    ddm2 = wtotal*wg2*ne2 * 1000
    ddma = wtotal*wga*nea * 1000
    ddmgc = wtotal*wggc*negc * 1000
    ddmlism = (1.-wgvN*wvoid)*wglism*wlism*nelism * 1000
    ddmcN = wgcN*necN * 1000
    ddmvN = wgvN*wvoid*nevN * 1000

    # cumulative integrals of n_e components:
    dm1run = array([trapz(ddm1[0:j], sf_vec[0:j]) for j in range(size(ddm1))])
    dm2run = array([trapz(ddm2[0:j], sf_vec[0:j]) for j in range(size(ddm2))])
    dmarun = array([trapz(ddma[0:j], sf_vec[0:j]) for j in range(size(ddma))])
    dmgcrun = array([trapz(ddmgc[0:j], sf_vec[0:j]) for j in range(size(ddmgc))])
    dmlismrun = array([trapz(ddmlism[0:j], sf_vec[0:j]) for j in range(size(ddmlism))])
    dmcNrun = array([trapz(ddmcN[0:j], sf_vec[0:j]) for j in range(size(ddmcN))])
    dmvNrun = array([trapz(ddmvN[0:j], sf_vec[0:j]) for j in range(size(ddmvN))])


    # spline functions for cumulative DM components:
    cs_dm1 = CubicSpline(sf_vec, dm1run)
    cs_dm2 = CubicSpline(sf_vec, dm2run)
    cs_dma = CubicSpline(sf_vec, dmarun)
    cs_dmgc = CubicSpline(sf_vec, dmgcrun)
    cs_dmlism = CubicSpline(sf_vec, dmlismrun)
    cs_dmcN = CubicSpline(sf_vec, dmcNrun)
    cs_dmvN = CubicSpline(sf_vec, dmvNrun)

    # DM component values at estimated distance
    dm1 = cs_dm1(dhat)
    dm2 = cs_dm2(dhat)
    dma = cs_dma(dhat)
    dmgc = cs_dmgc(dhat)
    dmlism = cs_dmlism(dhat)
    dmcN = cs_dmcN(dhat)
    dmvN = cs_dmvN(dhat)

    #-----------------------------
    # Individual components of SM:
        # [Note these need to be multiplied by
        # sm_factor = 1.836 to get SM in kpc m^{-20/3}]
    #-----------------------------

    dsm1 = wtotal*wg1*ne1**2*F1
    dsm2 = wtotal*wg2*ne2**2*F2
    dsma = wtotal*wga*nea**2*Fa

    # cumulative integrals of SM components (still divided by sm_factor)
    sm1run = array([trapz(dsm1[0:j], sf_vec[0:j]) for j in range(size(dsm1))])
    sm2run = array([trapz(dsm2[0:j], sf_vec[0:j]) for j in range(size(dsm2))])
    smarun = array([trapz(dsma[0:j], sf_vec[0:j]) for j in range(size(dsma))])
    smgcrun = array([trapz(dsmgc[0:j], sf_vec[0:j]) for j in range(size(dsmgc))])
    smlismrun = array([trapz(dsmlism[0:j], sf_vec[0:j]) for j in range(size(dsmlism))])
    smcNrun = array([trapz(dsmcN[0:j], sf_vec[0:j]) for j in range(size(dsmcN))])
    smvNrun = array([trapz(dsmvN[0:j], sf_vec[0:j]) for j in range(size(dsmvN))])

    # spline functions for cumulative SM components:
    cs_sm1 = CubicSpline(sf_vec, sm1run)
    cs_sm2 = CubicSpline(sf_vec, sm2run)
    cs_sma = CubicSpline(sf_vec, smarun)
    cs_smgc = CubicSpline(sf_vec, smgcrun)
    cs_smlism = CubicSpline(sf_vec, smlismrun)
    cs_smcN = CubicSpline(sf_vec, smcNrun)
    cs_smvN = CubicSpline(sf_vec, smvNrun)

    # SM components at estimated distance (now with sm_factor):
    sm1 = cs_sm1(dhat) * sm_factor
    sm2 = cs_sm2(dhat) * sm_factor
    sma = cs_sma(dhat) * sm_factor
    smgc = cs_smgc(dhat) * sm_factor
    smlism = cs_smlism(dhat) * sm_factor
    smcN = cs_smcN(dhat) * sm_factor
    smvN = cs_smvN(dhat) * sm_factor

    #---------------
    # LISM analysis:
    #---------------
    # Calculate path lengths through LISM components:
    # Take into account the weighting hierarchy, LHB:LOOPI:LSB:LDR

    ldr_path = lhb_path = lsb_path = loopI_path = 0.
    ldr_dist = lhb_dist = lsb_dist = loopI_dist = 0.

    lhb_path = trapz(wlhb, sf_vec)
    loopI_path = trapz((1-wlhb) * wloopI, sf_vec)
    lsb_path = trapz((1-wlhb)*(1-wloopI) * wlsb, sf_vec)
    ldr_path = trapz((1-wlhb)*(1-wloopI)*(1-wlsb) * wldr, sf_vec)

    if lhb_path > 0:
        lhb_dist = trapz(wlhb*sf_vec, sf_vec) / lhb_path 
    if loopI_path > 0:
        loopI_dist = trapz((1-wlhb)*wloopI*sf_vec, sf_vec) / loopI_path 
    if lsb_path > 0:
        lsb_dist = trapz((1-wlhb)*(1-wloopI)*wlsb*sf_vec, sf_vec) / lsb_path
    if ldr_path > 0:
        ldr_dist = trapz((1-wlhb)*(1-wloopI)*(1-wlsb)*wldr*sf_vec, sf_vec) / ldr_path

    # --------------------
    # Spiral arm analysis:
    # --------------------
    # Calculate path lengths through spiral arms:
    """
    pathlengths: 
        whicharm = 0,5 (currently).
            1,4 for the equivalent of the TC93 arms
            5 for the local arm
            0 means interarm path segments
    note narmsmax = number of arms
    armpaths has one additional entry to include interarm path length, so
    narmsmax1 = narmsmax + 1 (set in config_ne2001p.py)
    """
    armpaths = zeros(narmsmax1)
    armdistances = zeros(narmsmax1)

    # Using for loop for greater code clarity 
    # (could use list comprehensions but expressions are clunky)
    for j in range(narmsmax1):
        inds = where(whicharm==j) 
        sinds = sf_vec[inds]
        winds = whicharm[inds] / max(j, 1)              # avoid div by zero for j=0
        armpaths[j] = trapz(winds, sinds)
        armdistances[j] = trapz(sinds*winds, sinds) / where(armpaths[j] > 0, armpaths[j], 1) 

    # Printing to f24 file:
    print('LISM path lengths (kpc) with weighting hierarchy LHB:LOOPI:LSB:LDR',
    file = f24)
    print('{}'.format(18*' '+ '  LHB     LoopI     LSB      LDR'.expandtabs(15)),
    file = f24)
    print(3*' '+  '{0:15s}{1:6.3f}{2:9.3f}{3:9.3f}{4:9.3f}'.format('Length', lhb_path, loopI_path, lsb_path, ldr_path), file=f24)
    print(3*' '+  '{0:15s}{1:6.3f}{2:9.3f}{3:9.3f}{4:9.3f}'.format('Mean Dist.', lhb_dist, loopI_dist, lsb_dist, ldr_dist), file=f24)

    print('\nFractional contributions to DM:', file=f24)
    print('  outer   inner    arms     gc    lism    clumps  voids       DM', file= f24)
    print('{:7.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:11.3f}'.format(
        dm1/dm_target, dm2/dm_target, dma/dm_target, dmgc/dm_target,
        dmlism/dm_target, dmcN/dm_target, dmvN/dm_target, dm_target), file=f24)

    print('\nFractional contributions to SM:', file=f24)
    print('  outer   inner    arms     gc    lism    clumps  voids       SM', file= f24)
    print('{:7.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:11.3e}'.format(
        sm1/sm, sm2/sm, sma/sm, smgc/sm,
        smlism/sm, smcN/sm, smvN/sm, sm), file=f24)

    print('\nPath lengths through spiral arms:', file=f24)
    print('  Arm      Mean Distance       Path Length    (arm=0 => interarm)', file=f24)

    for j in range(narmsmax1):
         print('  {:2d}{:18.3f}{:18.3f}'.format(j, armdistances[j], armpaths[j]), file=f24)

    f24.close()
    return sf_vec, ne, ne1, ne2, nea

# ----------------------------------------------------------------------

def plot_dm_along_LoS(
        dm_target, dhat, sf_vec, dm_cumulate_vec, relevant_clump_indices, dc, 
        which = 'dm2d', plot_dm_target=False, saveplot=True):

        """
        Plots DM(s) along line of sight

        which = 'dm2d' or 'd2dm': controls labeling
        """ 

        # Plot cumulative DM 
        fig = figure()
        ax = fig.add_subplot(111)
        subplots_adjust(bottom=0.15)
        plot(sf_vec, dm_cumulate_vec, label=r'$\rm Model \ DM$')

        # Plot clump DMs if any affect the LoS
        if np.size(relevant_clump_indices) > 0:
            clump_distances = dc[relevant_clump_indices]
            clump_dm_hats = interp(clump_distances, sf_vec, dm_cumulate_vec) 
            plot(clump_distances, clump_dm_hats, 'ro', 
                label=r'$\rm  DM\,@\ clump \ positions)$')
        axis(xmin=-0.05*sf_vec[-1])
        axis(ymin=-0.05*dm_cumulate_vec[-1])
        plot((axis()[0], dhat), (dm_target, dm_target), 'k--', lw=1)
        plot((dhat, dhat), (axis()[2], dm_target), 'k--', lw=1)
        xlabel(r'$\rm Distance \ along \ LoS \ \ \  (kpc)$', fontsize=13)
        ylabel(r'$\rm DM(s) \ \    (pc\ cm^{-3})$', fontsize=13)

        if which == 'dm2d':
            title(r'$\rm DM\,\to\, D: \ \ DM_{target} = %5.1f \ pc\ cm^{-3} \quad\quad\quad \widehat d = %5.2f \ kpc$'%(dm_target, dhat), fontsize=12)
        if which == 'd2dm':
            title(r'$\rm D\,\to\, DM: \ \ D_{target} = %6.2f \ kpc  \quad\quad\quad \widehat DM(d_{\rm target})  = %5.1f \ pc\,cm^{-3}$'%(dhat, dm_cumulate_vec[-1]), fontsize=12)
            if plot_dm_target:
                plot(dhat, dm_target, 'ko', label=r'$\rm Target \ DM$')
        legend(loc=0, fontsize=10)
        annotate(plotstamp, xy=(0.70, 0.02), xycoords='figure fraction', ha='left', va='center', fontsize=5)

        savedir = os.getcwd()

        if saveplot:
            #plotfile = 'dm_vs_d_' + basename + '.pdf'
            plotfile = savedir+'/output_ne2001p/dm_vs_d_NE2001p.pdf'
            savefig(plotfile)
        show()
        
        return

# ----------------------------------------------------------------------

def plot_dm_sm_along_LoS(
    dm_target, dm_cumulate_vec, dhat, sf_vec,
    sm_hat, smtau_hat, smtheta_hat, smiso_hat,
    dsm, sm_cumulate, smtau_cumulate, smtheta_cumulate, smiso_cumulate,
    which = 'dm2d', saveplot=True):
    """
    Plots DM(s) and SM(s) along line of sight
    """

    # Plot cumulative DM and SMs
    fig = figure(figsize=(6,6))
    ax = fig.add_subplot(311)
    plot(sf_vec, sm_factor * dsm)
    ylabel(r'$\rm \Delta SM$',fontsize=13)
    if which == 'dm2d':
        title(r'$\rm DM\to D: \ \ \hat d = %5.2f \ \ SM, SM_\tau, SM_\theta = %8.4f \  %8.4f \ %8.4f$' %(dhat, sm_hat, smtau_hat, smtheta_hat), fontsize=10) 
    if which == 'd2ddm':
        title(r'$\rm D\to DM: \ \ \hat d = %5.2f \ \ SM, SM_\tau, SM_\theta = %8.4f \  %8.4f \ %8.4f$' %(dhat, sm_hat, smtau_hat, smtheta_hat), fontsize=10) 

    ax = fig.add_subplot(312)
    plot(sf_vec, sm_cumulate, label='SM')
    plot(sf_vec, smtau_cumulate, label=r'SM$_{\tau}$')
    plot(sf_vec, smtheta_cumulate, label=r'SM$_{\theta}$')
    plot(sf_vec, smiso_cumulate, label=r'SM$_{\rm iso}$')

    plot(dhat, sm_hat, 'o')
    plot(dhat, smtau_hat, 'o')
    plot(dhat, smtheta_hat, 'o')
    plot(dhat, smiso_hat, 'o')
    legend(loc=0)
    ylabel(r'$\rm SM(s)$',fontsize=13)

    ax = fig.add_subplot(313)
    plot(sf_vec, dm_cumulate_vec)
    xlabel(r'$\rm Distance \ along \ LoS \ \   (kpc)$', fontsize=13)
    ylabel(r'$\rm DM(s) \ \    (pc\ cm^{-3})$', fontsize=13)
    annotate(plotstamp, xy=(0.70, 0.02), xycoords='figure fraction', ha='left', va='center', fontsize=5)
    #show()

    fig.tight_layout()

    savedir = os.getcwd()
    
    if saveplot:
        plotfile = 'sm_and_dm_d_NE2001p.pdf'
        savefig(savedir+'/output_ne2001p/'+plotfile)
    show()

# ----------------------------------------------------------------------

# Main

if __name__ == '__main__':
    
    ldeg, bdeg, dm_target = 30, 0, 1000
    ndir = -1
    ds_fine = 0.005 # 0.01
    ds_coarse = 0.1

    l = deg2rad(ldeg)
    b = deg2rad(bdeg)

    verbose = True
    do_analysis = True
    dm2d_only = True
    dm2d_only = False
    plotting = False
    plotting = True
    debug = False 


    if dm2d_only:
        limit, dhat, dm_target \
        = dmdsm_dm2d(l, b, dm_target, ds_fine=ds_fine, ds_coarse=ds_coarse, 
              Nsmin=10, dm2d_only=dm2d_only, do_analysis=do_analysis, 
              plotting=plotting, verbose=verbose, debug=debug)
    else:
        limit, dhat, dm_target, sm, smtau, smtheta, smiso \
        = dmdsm_dm2d(l, b, dm_target, ds_fine=ds_fine, ds_coarse=ds_coarse, 
            Nsmin=10, dm2d_only=dm2d_only, do_analysis=do_analysis, 
            plotting=plotting, verbose=verbose, debug=debug)

    if plotting:
        input('hit return')
        close('all')

