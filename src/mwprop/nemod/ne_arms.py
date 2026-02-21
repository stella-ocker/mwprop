# mwprop.nemod v2.0 Jan 2026

"""
ne_arms
Replacement for Fortran routine NE_ARMS_LOG_MOD in density.NE2001.f

Version history:
2020 Jan 19
01/19/20 --- JMC
     * initial conversion f77 --> python
02/08/20 --- JMC
     * removed rsun = 8.5
     * all model parameters now imported from config_ne2001p
"""

from mwprop.nemod.config_nemod import *

script_path = os.path.dirname(os.path.realpath(__file__))

def ne_arms_ne2001p(x,y,z, Ncoarse=20, dthfine=0.01, nfinespline=5, verbose=False):
    """
    Evaluates electron density and F parameter for arm nearest to x,y,z

    Input:
        x, y, z = location in NE2001  coordinates  in kpc
        Ncoarse = number of coarse samples along each arm
        dthfine = step size in Galactocentric angle for fine sampling (rad)
        nfinespline = number of coarse samples to use for fine-sampling spline
        verbose = True:  writes out n_e component information

    Output:
        ne        electron density        cm^{-3}
        F         F parameter             composite units
        whicharm  which spiral arms       1 to 5, no close arm ==> 0
    """

    # First find coarse location of points on arms nearest to input point

    narms = np.shape(coarse_arms)[1]
    dsq_coarse = (coarse_arms[0] - x)**2 + (coarse_arms[1] - y)**2
    index_dsqmin = np.argmin(dsq_coarse, axis=1)

    # Zoom in to find precision location and evaluate density and F parameter

    # number of coarse samples to use as input to fine spline
    nspline2 = int((nfinespline-1)/2)

    # Initialize
    nea = 0.
    ga = 0.
    whicharm = 0
    whicharm_spiralmodel = 0

    # Cylindrical radius and Galactocentric angle of input point;
    # thxydeg measured CCW from +y axis

    rr = sqrt(x**2 + y**2)
    thxydeg = np.rad2deg(np.arctan2(-x, y))
    if thxydeg < 0:
        thxydeg += 360

    # Cache model parameters to avoid repeated dict lookups in the loop
    Dgal_wa = Dgal['wa']
    Dgal_Aa = Dgal['Aa']
    Dgal_ha = Dgal['ha']
    Dgal_Fa = Dgal['Fa']

    arm_index = np.fromiter((Darmmap[str(j)] for j in range(narms)), dtype=int)
    warm = np.array([Dgal['warm'+str(jj)] for jj in arm_index])
    harm = np.array([Dgal['harm'+str(jj)] for jj in arm_index])
    narm = np.array([Dgal['narm'+str(jj)] for jj in arm_index])

    if verbose:
        print('arms rr, thxydeg: ', rr, thxydeg)

    for j in range(narms):

        # Use fine spline on coarse samples  to find nearest position on j-th arm

        ind1 = range(max(0, index_dsqmin[j]-nspline2-1),
                     min(index_dsqmin[j]+nspline2+1, Ncoarse), 1)
        dsqs = CubicSpline(th1[j, ind1], dsq_coarse[j, ind1])
        thjfine = np.arange(th1[j, ind1[0]], th1[j, ind1[-1]], dthfine)
        ind_min = dsqs(thjfine).argmin()
        thjmin = thjfine[ind_min]
        # Use precomputed spiral arm spline when available to avoid per-call setup cost
        if 'armsplines' in globals() and len(armsplines) > j:
            sj = armsplines[j]
        else:
            sj = CubicSpline(th1[j,:], r1[j,:])
        rjmin = sj(thjmin)
        xjmin, yjmin = -rjmin*np.sin(thjmin), rjmin*np.cos(thjmin)

        # Evaluate electron density for nearest spiral arm if it is within
        # some multiple of the e-folding half-width of the arm

        dmin = sqrt((x-xjmin)**2 + (y-yjmin)**2)
        jj = arm_index[j]
        wa = Dgal_wa

        """
        Note armmap used here is to maintain the legacy arm numbering
        from NE2001 and TC93 using as input the arm numbering in Wainscoat.
        New NE model could dispense with this.
        """

        if j== 0: dmin_min = dmin
        Wa = wa * warm[j]
        if dmin < 3.*wa:
        #if dmin < 5.*wa:               # This helps reduce discontinuities
            if dmin <= dmin_min:
                dmin_min = dmin
                whicharm_spiralmodel = j+1

            # arm-distance factor:
            argxy = dmin / Wa
            ga = np.exp(-argxy**2)

            # Galactocentric radial factor:
            if rr > Dgal_Aa: ga *= sech2((rr-Dgal_Aa)/2.)

            # z factor:
            Ha = Dgal_ha * harm[j]
            ga *= sech2(z/Ha)

            # Amplitude re-scalings as in NE2001 code;
            # Ignore cases where these have been turned off in that code (fossils!).

            # TC arm 3:
            if jj == 3:
                th3adeg, th3bdeg = 290, 363
                test3 = thxydeg - th3adeg
                if test3 < 0: test3 += 360
                if 0 <= test3 and test3 < th3bdeg-th3adeg:
                    arg = 2.*pi*test3 / (th3bdeg-th3adeg)
                    fac = ((1 + np.cos(arg))/2)**4
                    # TEMP: fac = 1
                    ga *= fac

            # TC arm 2:
            if jj == 2:
                th2adeg, th2bdeg = 340, 370
                test2 = thxydeg - th2adeg
                fac2min = 0.1
                if test2 < 0: test2 += 360
                if 0 <= test2 and test2 < th2bdeg-th2adeg:
                    arg = 2.*pi*test2 / (th2bdeg-th2adeg)
                    fac = ((1 + fac2min + (1 - fac2min) * np.cos(arg))/2)**3.5
                    # fac = 1    # TEMP
                    ga *= fac

            nea += ga * narm[j] * Dgal['na']
            s = sqrt(x**2 + (rsun-y)**2)
            #s = sqrt(x**2 + (8.5-y)**2)
            if verbose:
                print('arms: ', j, jj, whicharm_spiralmodel, index_dsqmin[j],
                    max(0, index_dsqmin[j]-nspline2-1),
                    min(index_dsqmin[j]+nspline2+1, Ncoarse), ind1)
                print(
                    '%5.2f  %5.2f  %d  %d  %d    %4.2f  %4.2f  %4.2f    %4.2f  %7.5f'
                     %(s, rr, j, jj, whicharm_spiralmodel, wa, Wa, dmin, ga, nea))
    if whicharm_spiralmodel == 0:
        whicharm = 0
        Farm = 0
    else:
        whicharm = Darmmap[str(whicharm_spiralmodel-1)]
        #Farm = Dgal['Fa'] * Dgal['farm'+str(whicharm)]    # Intended value

        ###### NOTE ######
        # 2020 Feb 9:
        # found that the Fortran code doesn't use Farm calculated using
        # Fa x farm_j parameters.  This seems to be an error in dmdsm.NE2001.f,
        # which uses Fa to calculate SM for the spiral arms instead of Faval,
        # the value returned by density.NE2001.f.    This only affects the
        # spiral arm values for SM.

        # For now, maintain this error in the Python code because the aim here
        # is to replicate the Fortran code:

        Farm = Dgal_Fa

        # The question is then whether the error was in the code when fitting
        # was done to find the best values of farm_j?   I think the error was
        # not there during the fitting because some of the earlier code versions
        # use Fa * farm_j.

    if verbose:
        print('arms whicharm_spiralmodel  whicharm: ',
            whicharm_spiralmodel, whicharm)

    return nea, Farm, whicharm