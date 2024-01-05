# mwprop.ne2001p v1.0 Jan 2024

'''
Python version of subroutine nevoidN.f in NE2001 Fortran code 

Returns electron density nevN and fluctuation parameter FvN 
at position designated by l,b,d,x,y,z c for a set of  
voids with parameters read in from file  nevoidN.dat

input:
    x,y,z   coordinates (kpc)  (as in TC93)

output:
    nevN    electron density in void at (x,y,z)
    FvN fluctuation parameter
    hitvoid =   0:     no void hit
                j>0:   j-th void hit
    wvoid = 0,1:     void weight    

parameters:
    lv  = galactic longitude of void center
    bv  = galactic latitude of void center
    dv  = distance from Sun of void center
    (xv,yv,zv) = void center location (calculated)
        nev = internal peak electron density
        Fv      = void fluctuation parameter
    aav = void major axis at 1/e
    bbv = void minor axis at 1/e
    ccv = void minor axis at 1/e
    thvy    = rotation axis of void about y axis
    thvz    = rotation axis of void about z axis
    edgev   = 0 => use exponential rolloff out to 5rc
                 1 => uniform and truncated at 1/e
Version history:

01/20/20 Stella Koch Ocker 
    * initial conversion f77 --> python
01/23/20 -- now reads input parameters from dictionary program ne2001p_input
02/08/20 -- JMC
    * imports config_ne2001p for  model setup
    * rsun = 8.5 commented out; now set in setup file
02/09/20 -- JMC
    * comment out hitvoidflag; not used even in Fortran code
02/10/20 -- JMC
    * definition of void arrays moved to config_ne2001p.py
'''

from mwprop.ne2001p.config_ne2001p import *

def nevoidN(x,y,z):

    nevN = 0.
    FvN = 0.
    hitvoid = 0
    wvoid = 0

    '''
    note rotation matrix in the 'q = ' statement below
    corresponds to \Lambda_z\Lambda_y
    where \Lambda_y = rotation around y axis
        \Lambda_z = rotation around z axis
    defined as
        \Lambda_y =  c1  0  s1
                     0  1   0
                    -s1  0  c1

        \Lambda_z =  c2 s2   0
                    -s2 c2   0
                    0  0   1
        =>
        \Lambda_z\Lambda_y =  c1*c2   s2   s1*c2
                            -s2*c1   c2  -s1*s2
                            -s1    0      c1
    so the rotation is around the y axis first, then the z axis
    '''

    for j in range(nvoids):
        dx = x-xv[j]
        dy = y-yv[j]
        dz = z-zv[j]
        q = (cc12[j]*dx + s2[j]*dy + cs21[j]*dz)**2. / aav[j]**2. + (-cs12[j]*dx + c2[j]*dy - ss12[j]*dz)**2. / bbv[j]**2. + (-s1[j]*dx + c1[j]*dz)**2. / ccv[j]**2.
        if edgev[j] == 0. and q < 3.: # note this doesn't actually get used; no clumps with edge = 0
            nevN = nev[j] * exp(-q)
            FvN = Fv[j]
            hitvoid = j+1
            #hitvoidflag[j] = 1
        if edgev[j] == 1. and q <= 1.:
            nevN = nev[j]
            FvN = Fv[j]
            hitvoid = j+1 # 3/6/22 -- SKO changed this from hitvoid = j --- j = 0 for Gum edge, which means it doesn't get counted
            #hitvoidflag[j] = 1

    if hitvoid != 0:
        wvoid = 1
        
    return nevN, FvN, hitvoid, wvoid
