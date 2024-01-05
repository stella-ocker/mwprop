#!/usr/bin/env python

# mwprop.ne2001p v1.0 Jan 2024


"""
los_diagnostics

Based on earlier los_plot.py and other code
Can vary Ncoarse = number of coarse samples along the LoS to test precision:
    Ncoarse  = 20           # as in Fortran code, produces discontinuities
    Ncoarse = 100           # reduces discontinuities, but no better than Ncoarse=50
    Ncoarse = 50            # reduces discontinuities

JMC 2024 Jan 3 
"""

from matplotlib.pyplot import suptitle

import mwprop.ne2001p.NE2001 as NE2001

import mwprop.ne2001p.config_ne2001p as config_ne2001p
import argparse

from mwprop.ne2001p.density_components import *
import mwprop.ne2001p.ne_arms_ne2001p as ne_arms 

import mwprop.ne2001p.density_components as dc
import mwprop.ne2001p.density_ne2001 as ne01
import mwprop.ne2001p.ne_lism as lism
import mwprop.ne2001p.neclumpN_NE2001_fast as necf
import mwprop.ne2001p.nevoidN_NE2001 as nev

script_path = os.path.dirname(os.path.realpath(__file__))
basename = sys.argv[0].split('/')[-1].split('.')[0]         # for plot file name

def plot_ne_arms():
    # Plotting
    # ------------------------------------------------------------------
    # Three panel plot with n_e vs s, spiral arms, and arm number vs LoS
    # ------------------------------------------------------------------
    fig = figure()
    subplots_adjust(hspace=0.35, left=0.15, bottom=0.15, top = 0.85)
    suptitle(r'$\rm Output \ from \ los\_diagnostics.py$', fontsize=12)

    # Frame 1: ne from spiral arms vs distance
    ax1 = fig.add_subplot(211)
    plot(s25, ne25, 'k-', lw=6, label=r'$\rm n_e\, (total) $', alpha=0.7)
    plot(s, ne_lism_vec, 'r-', lw=3, label=r'$\rm {n_e}_{,lism} $')
    plot(s, nea_vec, '-', ms=2, label=r'$\rm n_e\, (arm) $')
    plot(s, ne1_vec, '-', ms=2, label=r'$\rm {n_e}_1 $')
    plot(s, ne2_vec, '-', ms=2, label=r'$\rm {n_e}_2 $')
    #plot(s, negc_vec, '-', ms=2, label=r'$\rm {n_e}_{gc} $')
    xlabel(r'$\rm d \ \ (kpc) $', fontsize=12)
    ylabel(r'$\rm n_e \ \ (cm^{-3}) $', fontsize=12)
    yy = axis()[3]
    axis(ymax = 1.35*yy)
    legend(loc=(0.025, 0.8), ncol=5, fontsize=9)
    if ndir > 0:
        title(r'$\rm DM = %5.1f\ pc\,cm^{-3} \ \longrightarrow \   D = %5.1f\ kpc \ \ \ \ \ \ \ \  l,\ b \ = \ %6.2f, \   %6.2f$'%(DM, d, ldeg, bdeg), fontsize=11)
    else:
        title(r'$\rm D = %5.1f\, kpc  \ \longrightarrow \   DM = %5.1f\, pc\ cm^{-3} \ \ \ \ \ \ \ \  l,\ b \ = \ %6.2f, \   %6.2f$'%(d, DM, ldeg, bdeg), fontsize=11)

    # Frame 2: spiral arms
    ax2 = fig.add_subplot(223)
    ax2.set_aspect('equal',  'box')
    for j in range(Darms['a'].size): 
        plot(armarray[0, j], armarray[1, j])
        #plot(coarse_arms[0, j], coarse_arms[1, j])
    # line of sight
    # Sun location
    plot(0, rsun, 'o', lw=0.5, ms=4, mfc='w', mec='k')
    plot(0, rsun, '.', ms=1, mfc='k', mec='k')

    # Galactic center
    plot(0, 0, '+')

    # line of sight
    plot(xvec, yvec, 'k-', lw=1)

    xlabel(r'$\rm X \ \ (kpc) $', fontsize=12)
    ylabel(r'$\rm Y \ \ (kpc) $', fontsize=12)

    # Frame 3: arm numbers
    ax3 = fig.add_subplot(224) 
    for j in range(1, Darms['a'].size+1):
       inds = np.where(Armvec==j)
       plot(s[inds], Armvec[inds], lw=2)
    xlabel(r'$\rm d \ \ (kpc) $', fontsize=12)
    ylabel(r'$\rm Arm \ number  $', fontsize=12)

    # save file in main
    #show()
    return
   
def plot_dm_ne_cn2():
    # ----------------------------------
    # Three panels with DM, n_e, and Cn2
    # ----------------------------------
    fig = figure()
    subplots_adjust(left=0.20, bottom=0.15)
    suptitle(r'$\rm Output \ from \ los\_diagnostics.py$', fontsize=12)

    # Frame 1:  DM vs s
    ax = fig.add_subplot(311)
    plot(s25, dm25)
    ylabel(r'$\rm DM \ (pc\ cm^{-3}) $', fontsize=14, labelpad=30)
    annotate(r'$\rm l = %6.1f^{\circ}  \ \ b = %5.1f$'%(ldeg, bdeg), xy=(0.025, 0.875), xycoords='axes fraction', ha='left', va='center', fontsize=10)
    annotate(r'$\rm DM_{max} =  %8.2f$'%(DMmax), xy=(0.025, 0.650), xycoords='axes fraction', ha='left', va='center', fontsize=10)
    title(r'$\rm Effective \ distances: \ \overline{d}_{n_e} = %8.2f \ \ \ \ \ \overline{d}_{n_e^2} = %8.2f \ \ \ \ \ \overline{d}_{SM} = %8.2f$'%(dbar_ne, dbar_ne2, deffsm2), fontsize=11)
    tick_params(labelbottom = False)

    # Frame 2: n_e vs s
    ax = fig.add_subplot(312)
    plot(s25, ne25, label=r'$\rm n_e$')
    ylabel(r'$\rm n_e \ (cm^{-3}) $', fontsize=14, labelpad=20)
    tick_params(labelbottom = False)

    # Frame 3: Cn2e vs s
    ax = fig.add_subplot(313)
    plot(s25, Cn2_vs_s25, label=r'$\rm C_n^2$')
    ylabel(r'$\rm  C_n^2 \ (m^{-20/3}) $', fontsize=14, labelpad=10)

    xlabel(r'$\rm Distance \ from \ solar \ system \ \ (kpc) $', fontsize=16)
    #legend(loc=0)
    #show()

    return

# Main

if __name__ == '__main__':

    try:
        parser = argparse.ArgumentParser(
            description='Plot NE2001 diagnostics for designated line of sight)')

        """
        parser.add_argument('-l', '--l', default=57.51, help='Galactic longitude (deg)')
        parser.add_argument('-b', '--b', default=-0.29, help='Galactic latitude (deg)')
        parser.add_argument('-DM_D', '--DM_D', default=71, help='DM (pc cm^{-3}) or distance (kpc)')
        parser.add_argument('-ndir', '--ndir', default=1, help='ndir = 1 (DM->D), ndir = -1 (D->DM)')
        """

        parser.add_argument('l', help='Galactic longitude (deg)')
        parser.add_argument('b', help='Galactic latitude (deg)')
        parser.add_argument('DM_D', help='DM (pc cm^{-3}) or distance (kpc)')
        parser.add_argument('ndir', help='ndir = 1 (DM->D), ndir = -1 (D->DM)')

        # Default values for sampling along the line of sight: less likely to need change. 
        Ncoarse_def = 50
        ds_def = 0.01
        parser.add_argument('-N', '--Ncoarse', default=Ncoarse_def, 
                help='Number of coarse values to use for smooth n_e components')
        parser.add_argument('-ds', '--ds', default=ds_def, 
                help='step size along line of sight (pc)')

        args = parser.parse_args()

        ldeg = float(args.l)
        bdeg = float(args.b)
        dmd = float(args.DM_D)
        ndir = int(args.ndir)
        Ncoarse = int(args.Ncoarse)
        ds = float(args.ds)

        l, b = deg2rad((ldeg, bdeg))

        # Directory for reading the f25 file written by NE2001 (from call below)
        f25dir = os.getcwd()+'/output_ne2001p/'

        # Plot directory (same as f25dir here):
        plotdir = os.getcwd()+'/output_ne2001p/'

        print('Line of sight diagnostics for l = %6.1f   b = %6.1f   dmd = %d   ndir = %d'%(ldeg, bdeg, dmd, ndir))
        if Ncoarse != Ncoarse_def:
            print('Using Ncoarse = %d (number of samples used for smooth n_e components)'%(Ncoarse)) 
        if ds != ds_def:
            print('Using ds = %6.3f (fine sample interval for LoS)'%(ds)) 
        print('Plots are stored as PDFs in directory: \n%s'%(plotdir))

    except SystemExit:
        print('Try again with inputs:')
        print('    Use  los_diagnostics.py  l       b     DM_D ndir   with DM_D = DM or distance for ndir > or < 0')
        print('    e.g. los_diagnostics.py 57.51 -0.29   71.02   1    for B1937+21 parameters')
        sys.exit()

    # Get output DM or d  from NE2001:
    xx = NE2001.ne2001(ldeg, bdeg, dmd, ndir, classic=False, dmd_only=False, do_analysis=True)

    if ndir > 0:                # DM -> D
        DM = xx[1]['DM']
        D2001 = xx[1]['DIST']
        d = D2001
    else:
        DM2001 = xx[1]['DM']
        d = xx[1]['DIST']
        DM = DM2001

    # get spiral arm info for plotting
    Dgal01, Dgc, Dlism, Dclumps, Dvoids, Darms, Darmmap, armmap, \
           r1, th1, th1deg, coarse_arms, rarray, tharray,  armsplines, armarray, \
           tangents, normals, curvatures = config_ne2001p.setup_spiral_arms(Ncoarse=Ncoarse)

    ne_arms.th1 = th1
    ne_arms.r1 = r1
    ne_arms.coarse_arms = coarse_arms
    ne_arms.Darmmap = Darmmap
    ne_arms.armmap = armmap
    ne_arms.Dgal01 = Dgal01
    
    s = np.linspace(0, d, int(d/ds))
    xvec, yvec, zvec =  s*np.sin(l), rsun-s*np.cos(l), s*np.sin(b)  

    nea_vec = zeros(np.size(xvec))
    ne1_vec = zeros(np.size(xvec))
    ne2_vec = zeros(np.size(xvec))
    negc_vec = zeros(np.size(xvec))
    Fvec = zeros(np.size(xvec))
    Armvec = zeros(np.size(xvec), dtype=int)
 
    ne_lism_vec = zeros(np.size(xvec))
    F_lism_vec = zeros(np.size(xvec))
    wlism_vec = zeros(np.size(xvec))

    for n, x in enumerate(xvec):

        # Get individual electron density contributions
        nea_vec[n], Fvec[n], Armvec[n] = \
            ne_arms.ne_arms_ne2001p(x ,yvec[n], zvec[n], Ncoarse=Ncoarse, verbose=False)
        ne1_vec[n] = dc.ne_outer(x, yvec[n], zvec[n])[0] 
        ne2_vec[n] = dc.ne_inner(x, yvec[n], zvec[n])[0] 
        negc_vec[n] = dc.ne_gc(x, yvec[n], zvec[n])[0] 

        ne_LISM, FLISM, wLISM, wLDR, wLHB, wLSB, wLOOPI = lism.ne_lism(xvec[n], yvec[n], zvec[n])
        ne_lism_vec[n] = ne_LISM
        F_lism_vec[n] = FLISM
        wlism_vec[n] = wLISM

    # Get total n_e from file f25 text file written by analysis_dmd_dm_and_sm() in dmdsm_ne2001p.py
    if ndir < 0:
        f25file = 'f25_d2dm_ne_dsm_vs_s.txt'
    else:
        f25file = 'f25_dm2d_ne_dsm_vs_s.txt'
        
    # Read f25 output file written by NE2001 into standard directory:
    ldeg, bdeg = np.loadtxt(f25dir+f25file, unpack=True, usecols = (0, 1), skiprows=1, max_rows=1)
    s25, ne25, Cn2_vs_s25, dm25, nea25 = np.loadtxt(f25dir + f25file, unpack=True, usecols = (0, 4, 5, 10, 11), skiprows=3)

    DMmax = dm25.max()

    # Effective distances based on LoS integrals weighted by various functions of n_e:
    dbar_ne = np.trapz(s25*ne25, s25) / np.trapz(ne25, s25)
    dbar_ne2 = np.trapz(s25*ne25**2, s25) / np.trapz(ne25**2, s25)
    deffsm2 =  np.trapz(s25*Cn2_vs_s25, s25) / np.trapz(Cn2_vs_s25, s25)

    # Cumulative DM and SM: not currently used, so commented out; may use later so keep.
    #DMvec = 1000.*np.array([np.trapz(ne25[:j], s25[:j]) for j in range(np.size(s25))])
    #SMvec = np.array([np.trapz(Cn2_vs_s25[:j], s25[:j]) for j in range(np.size(s25))])


    plot_ne_arms()
    plotfile = 'ne_vs_d_arms_' + '%6.1f'%(ldeg) + '_' + '%5.1f'%(bdeg) + '_' + str(int(DMmax)) + '_' +  basename + '.pdf'
    plotfile = plotfile.replace(' ', '')
    savefig(plotdir+plotfile)
    show()

    plot_dm_ne_cn2()
    plotfile = 'dm_ne_cn2_l_b_dmmax_' + '%6.1f'%(ldeg) + '_' + '%5.1f'%(bdeg) + '_' + str(int(DMmax)) + '_' +  basename + '.pdf'
    plotfile = plotfile.replace(' ', '')
    savefig(plotdir+plotfile)
    show()
    
    input('hit return')
    close('all')
