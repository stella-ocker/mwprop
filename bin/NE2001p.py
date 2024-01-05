#!/usr/bin/env python

# mwprop.ne2001p v1.0 Jan 2024

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from mwprop import ne2001p
from mwprop.ne2001p.NE2001 import *


if __name__ == '__main__':

    # If 'explain' option is set, print out README file and exit.
    if '-e' in sys.argv or '--explain' in sys.argv:
        try:
            inp_file = (pkg_resources.files(ne2001p) / 'README_NE2001p.txt')
            with inp_file.open('rt') as f:
                template = f.read()
                print(template)
        except AttributeError:
            template = pkg_resources.read_text(ne2001p, 'README_NE2001p.txt')
            print(template)
        sys.exit()

    try:
        # parse command line arguments and execute
        parser = argparse.ArgumentParser(description='NE2001p v1.0 Jan 2024')
        parser.add_argument('l',help='Galactic longitude (deg)')
        parser.add_argument('b',help='Galactic latitude (deg)')
        parser.add_argument('DM_D',help='DM (pc cm^{-3}) or distance (kpc)')
        parser.add_argument('ndir',help='ndir = 1 (DM->D), ndir = -1 (D->DM)')

        parser.add_argument('-a', '--analysis', action='store_true', 
            help='Do line of sight analysis')

        parser.add_argument('-p', '--plotting', action='store_true', 
            help='Do diagnostic plotting')

        parser.add_argument('-s', '--scattering', action='store_true', 
            help='Calculate scattering and scintillation parameters')

        parser.add_argument('-m', '--modern', action='store_true', 
            help='Modern output (turns off classic output like Fortran version)')

        parser.add_argument('-v', '--verbose', action='store_true', 
            help='Verbose output (classic output of Fortran version)')

        parser.add_argument('-b', '--debug', action='store_true', 
            help='debug output')
            
        args = parser.parse_args()


    except SystemExit:
        print('Try again with inputs')
        print('Use NE2001p.py -e to get explanation of code')
        sys.exit()

    ldeg = float(args.l)
    bdeg = float(args.b)
    dmd = float(args.DM_D)
    ndir = int(args.ndir)
    do_analysis = args.analysis
    do_plotting = args.plotting
    calc_scattering = args.scattering
    dmd_only = not calc_scattering
    verbose = args.verbose
    debug = args.debug
    modern = args.modern
    
    if modern:
        classic = False
    else:
        classic = True

    if do_plotting and calc_scattering:
        do_analysis = True # required for scattering plots

    Dn, Dv, Du, Dd = ne2001(ldeg, bdeg, dmd, ndir, 
        classic=classic, verbose=verbose, dmd_only=dmd_only, 
        do_analysis=do_analysis, plotting=do_plotting, debug=debug)

    if calc_scattering and do_plotting:
        output_dir = os.getcwd()+'/output_ne2001p/'
        if ndir >= 0:
            plot_dm_ne_cn2_vs_d(Dv, f25file = output_dir+'f25_dm2d_ne_dsm_vs_s.txt')
        if ndir < 0: # SKO added 12-15-23 to make output plotting work for both DM->D and D->DM
            plot_dm_ne_cn2_vs_d(Dv, f25file = output_dir+'f25_d2dm_ne_dsm_vs_s.txt')

        input('hit return')
        close('all')
