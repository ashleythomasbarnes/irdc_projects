import numpy as np
import math
from astropy.io import fits
import astropy.units as units
import astropy.units as au
import astropy.constants as constant
from reproject import reproject_interp
import matplotlib.pyplot as plt

# define continuum to column density function
hplanck          = constant.h.value
clight           = constant.c.value
kboltzmann       = constant.k_B.value
pc2cm            = units.pc.to('cm')
msun             = units.Msun.to('g')
mh               = constant.u.value
mh_g             = constant.u.to('g').value
G                = constant.G.value
sin1yr           = units.yr.to('s')
pc2m             = units.pc.to('m')
percm2perm       = (units.m**3).to('cm^3')
# mh               = constant.u.to('g')

def planck_wave( Wave, Temp ):

    "Planck function using wavelength"

    planck_conv_wave = 1.e-26 * clight / Wave**2.0 # Conv to useful astronomical units

    planck = ((2.0*hplanck*clight**2.0)/(Wave**5.0))*(1.0/(np.exp((hplanck*clight)/(Wave*kboltzmann*Temp))-1.0))
    planck = planck/planck_conv_wave

    return planck

def remove_nans(arr1, arr2):
    "Module to remove nans from several arrays -- returns flattend array."

    arr1 = np.asarray( arr1.flatten() )
    arr2 = np.asarray( arr2.flatten() )

    nan_ID1 = np.isnan(arr1)
    arr1_nonnan = arr1[~nan_ID1]
    arr2_nonnan = arr2[~nan_ID1]

    nan_ID2 = np.isnan(arr2_nonnan)
    arr1_nonnan = arr1_nonnan[~nan_ID2]
    arr2_nonnan = arr2_nonnan[~nan_ID2]

    return arr1_nonnan, arr2_nonnan

def column_density(flux, beam, frequency, temp, beta = 1.75, gas2dust = 100, mu = 2.8):

    "General column density calculation, beam in arcsec"

    beam = np.radians(beam / 3600.)

    wavelength = constant.c.value / frequency

    kappa = (0.9 * (frequency / 230.e9) **beta)

    B = planck_wave( wavelength, temp )

    omega = (math.pi / (4.0 * math.log(2.)) ) * beam[0]*beam[1]

    N_h2 = (flux * gas2dust) / (omega * mu * (mh * 1.e3) * kappa * B)

    return N_h2


def mass_calc(flux, frequency, temp, dist, beta = 1.75, gas2dust = 100, mu = 2.8):  #Wave, Temp, Kappa, Integrated_Flux, Obj_Dist ):

    # A more general mass calculation

    wavelength = constant.c.value / frequency
    dist_cm = dist * pc2cm

    kappa = (0.899 * (frequency / 230.e9) **beta)

    B = planck_wave( wavelength, temp )

    Mass = (dist_cm**2. * flux * gas2dust) / (kappa * B)
    Mass = Mass / msun

    return Mass


def number_density_sphere_pc( Mass_sol, Radius_pc, mu = 2.8):

    Mass = Mass_sol * msun
    Radius = Radius_pc * pc2cm
    V = (4. / 3.) * np.pi * (Radius**3.0)
    M = Mass / (mu * mh_g)

    n = M / V

    return n


def get_rho(mass_sol, radius_pc):
    """Get volume density"""
    mass = (mass_sol*au.Msun).to('kg')
    radius = (radius_pc*au.pc).to('m')
    
    vol = (4./3.)*np.pi*(radius**3.0)
    rho = mass/vol
    return(rho)

def tff_spherical(number_density, mu = 2.8):

    mass_density = mu * mh * number_density * percm2perm
    tff = np.sqrt( (3. * np.pi) / (32. * G * mass_density) )
    tff = tff / sin1yr

    return tff /1e6


def get_histo(data, weights = 1., n_bins = 75, threshold = 3e23, binlims = ''):

    data_1d = data.ravel()
    data_1d_noneg = data_1d[np.where(data_1d>=0)]
    mean = np.nanmean(data_1d_noneg)
    x_crit_col = threshold / mean

    if binlims == '':
        bins = 10. ** np.linspace(np.log10(np.nanmin(data_1d_noneg) - (np.nanmin(data_1d_noneg) * 0.1) ),
                                  np.log10(np.nanmax(data_1d_noneg) + (np.nanmax(data_1d_noneg) * 0.1)), n_bins)
    else:
        bins = 10. ** np.linspace(np.log10(binlims[0]), np.log10(binlims[1]), n_bins)

    histo = np.histogram(data_1d_noneg[~np.isnan(data_1d_noneg)], bins=bins)

    weights = 1. / np.nanmax(histo[0])
    histo = np.array([np.append(histo[0], 0) * weights, histo[1]])
    histo[0 , np.where(histo[0] == 0)] = 1e-10

    return x_crit_col, bins, histo, weights


def mass_from_column(cdense, area, distance, mu=2.8):

    """
    The mean column denisty, area in arcsec, distance in pc
    """

    distance = distance * pc2cm
    conv_distance = distance / 206265.

    area_pc = area*(conv_distance**2)

    mmp = mu*mh

    mass = cdense * area_pc * mmp
    mass_sol = mass / (msun/1e3)

    return mass_sol
