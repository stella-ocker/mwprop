# Instructions to read in calibration data for  NE2025, in Python:

### Distances ###

# See Section 2 of the paper for details on where the distances are drawn from. N = 171 distances, all with fractional uncertainty <25%. 

# load globular cluster distances
gcdat = load('GCdata.npz')
gc_names = gcdat['name']
gc_ngc = gcdat['ngc_name']
gc_gl = gcdat['gl']
gc_gb = gcdat['gb']
gc_dm = gcdat['dm']
gc_dmerr = gcdat['dmerr'] # standard error of the mean for all pulsars within each cluster
gc_dist = gcdat['dkpc']
gc_derr = gcdat['dkpc_err']
gc_names[gc_names=='none'] = gc_ngc[gc_names=='none'] # fills gaps in catalog names

# load parallax distances
pxdat = load('psr_dist_data.npz')
px_jnames = pxdat['jnames']
px_gl = pxdat['gl']
px_gb = pxdat['gb']
px_dm = pxdat['dm']
px_dkpc = pxdat['dkpc']
px_deplus = pxdat['deplus']
px_deminus = pxdat['deminus']
px_eplus = pxdat['px_eplus']
px = pxdat['px']

# gather all distances into single arrays
all_gl = concatenate((gc_gl,px_gl)) # all longitudes
reordered_gl = copy(all_gl)
reordered_gl[all_gl>180]-=360 # longitudes >180 set to negative values
all_gb = concatenate((gc_gb,px_gb)) # all latitudes
all_dm = concatenate((gc_dm,px_dm)) # all DMs
all_dkpc = concatenate((gc_dist,px_dkpc)) # all distances
all_deplus = concatenate((gc_derr,px_deplus)) # all positive distance errors
all_deminus = concatenate((gc_derr,px_deminus)) # all negative distance errors
all_names = concatenate((gc_names,px_jnames)) # all object names

### Scattering Measurements ###

# Measurements from Cordes et al. (2022), plus new FAST measurements from Jing et al. (2025). FAST measurements are only included if they have reasonable spectral indices (consistent with 4 - 4.4). If there are multiple scattering measurements for a given pulsar, they've been averaged together.

taudat = load('scattering_database_nov2025.npz')
mean_tau = taudat['tau'] # ms at 1 GHz
mean_dm = taudat['dm']
tau_names = taudat['psr_names']
tau_gl = taudat['gl']
tau_gb = taudat['gb']