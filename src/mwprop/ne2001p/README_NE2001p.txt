MWPROP

Jan 2024 v1.0

MWPROP provides NE2001p, a native Python implementation of the original Fortran code for NE2001. The package also contains a required `scattering_functions` module. NE2001p is accessible from the command line, similar to the Fortran code, or within Python scripts (see below). Future package releases will include additional modules for extensions of the NE2001 model. 

Please cite the NE2001p research note (Ocker & Cordes 2024) for use of this package.
The NE2001 model is described in detail in Cordes & Lazio (2002; 2003).

-----

Installation:

With pip:

>> pip install mwprop

On GitHub: https://github.com/stella-ocker/mwprop

Executable scripts `NE2001p.py` and `los_diagnostics.py` are automatically installed with pip and may be run as executables from the command line in any directory.

Dependencies
- python >= 3.6 (might work with python >= 3.0)
- numpy
- matplotlib
- scipy
- astropy
- mpmath

-----

Comparison to Fortran Version:

The input parameters and model components for the Milky Way are identical. 

Output is nearly the same but not identical because integrations are done slightly differently to speed up the Python version. 

Differences:

1.  Scattering computations are done as an option 
    to give maximum speed for DM to D or D to DM computations
2.  Output can be printed in console or returned in Python dictionaries, Dn, Dv, Du, Dd:
        Dn = names of variables
        Dv = values 
        Du = units
        Dd = descriptions
3.  Numerical integrations use numpy.trapz instead of step-by-step summing.
4.  Clumps are prefiltered for inclusion along a specified line of sight (LoS).
5.  Large-scale components (thin and thick disks, spiral arms) are sampled coarsely
    (0.1 kpc default sample interval) and cubic splines are used to resample onto
    a fine grid. 
6.  Components with small scale structure (local ISM, Galactic center, clumps, voids)
    are sampled directly on a fine grid.
7.  Different routines are used for execution to find distance from DM and for DM from distance.
    This is transparent to the user. Cases where the DM exceeds the maximum Galactic value will return a warning; distances predicted for these cases are unreliable (see Ocker & Cordes 2024).
8.  Config file:  code imports,  reading parameter files, and setting key values
    (like the Sun-GC distance = 8.5 kpc) are done by importing mwprop.ne2001p.config_ne2001p  
    Note the 8.5 kpc distance is needed because model parameters were determined using
    this value for the original NE2001 code. It cannot be changed by itself.

NE2001p is about 45 times slower than the Fortran implementation. 

Forward plan:
1. NE2001p = this version (p for Python)
2. NE2001x = next version (x for extended)
3. NE202(...)  = new version, complete redefinition of components and implementation.

----- 

Usage:

Command Line Usage: NE2001p.py ldeg bdeg dmd ndir [-options]

Required arguments:  ldeg bdeg dmd ndir     (Same as Fortran version)

    ldeg, bdeg = Galactic longitude and latitude in deg
    dmd = DM (pc/cc) or d (kpc)
    ndir >= 0   DM -> d
          < 0    d -> DM

Optional requirements and defaults:

    -a  --analysis   Calculates contributions from different model components and writes to file in user's working directory.
                     [False]
    -p  --plotting   Plots DM, n_e, scattering vs distance along LoS; saves pdf files in user's working directory.
                     [False]
    -s  --scattering Calculates scattering and scintillation parameters, etc.
                     [False] 
    -m  --modern     Modern output (suppresses classic Fortran output)
                     [False]
    -v  --verbose    Writes results to std out in the 'classic' form of the Fortran version.
                     [False]
    -b  --debug      Prints various diagnostics for clumps and integration. Only works for DM->D.
                     [False]
    -e  --explain    Show long-form explanation of code.
                     [False]

Script/iPython Usage:

The ne2001() function evaluates NE2001p and can be imported from the mwprop.ne2001p.NE2001 module. The function output consists of four dictionaries: 

    Dk    => Dictionary keys for output values
    Dv    => Numerical output values
    Du    => Units of output values
    Dd    => Extended output description 

Additional options for ne2001():

    classic = True      => print traditional Fortran output (default = True)
    dmd_only = True     => compute only DM or distance, no scattering quantities (default = False)
    do_analysis = True  => calculate and output diagnostic quantities for the line of sight (default = False)
    plotting = True     => plot DM vs distance along LoS, SM also if dmd_only = False (default = False)
    verbose = True      => print results to std out (default = False)

iPython Example:

    >>> from mwprop.ne2001p import *
    >>> from mwprop.ne2001p.NE2001 import ne2001
    >>> Dk,Dv,Du,Dd = ne2001(ldeg=45,bdeg=5,dmd=50,ndir=1,classic=False)
    >>> Dk['DIST']
    'DIST'
    >>> Dv['DIST']
    2.6366
    >>> Du['DIST']
    'kpc'
    >>> Dd['DIST']
    'ModelDistance'

In analysis and plotting modes, output files are saved to a folder called 'output_ne2001p' in the user's working directory. 

-----
    
Diagnostic code los_diagnostics.py

Plots electron density, DM, and C_n^2 along the line of sight designated by l, b, DM or d, and ndir = 1 or -1 (as with NE2001).
Also shows the line of sight projected onto the Galactic plane along with spiral arms used in NE2001.

This code can be run from any directory if `mwprop` is fully installed. Outputs are saved to a folder created in the user's working directory.

    >>> Usage:
    >>> los_diagnostics.py  l  b  dmd  ndir   
    >>> with l, b in deg, dmd = DM (pc/cc) or distance (kpc) for ndir > or < 0
    >>> [same scheme as input to NE2001]
