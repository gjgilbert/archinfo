import numpy as np
import warnings
import scipy.optimize as op
import .LMC
from   scipy import stats

pi = np.pi

MSME = 332948.6   # (M_sun/M_earth)
BIGG = 6.67e-11    # Newton's constant [SI units]
RSUN = 6.957e8     # solar radius [m]
MSUN = 1.988e30    # Solar mass [kg]


###

__all__ = ["dynamical_mass", "mass_partitioning", "monotonicity", \
           "characteristic_spacing", "gap_complexity", "flatness"]

# Functions containing relevant physics

def P_to_a(P, Mstar):
    """
    Convenience function to convert periods to semimajor axis from Kepler's Law
    
    Parameters
    ----------
    P : array-like
        orbital periods [days]
    Mstar : float
        stellar mass [solar masses]
        
    Returns
    -------
    a : array-like
        semi-major axis [stellar radii]
    """
    Pearth = 365.24    # [days]
    aearth = 215.05    # [solar radii]
    
    return aearth * ((P/Pearth)**2 *(1/Mstar))**(1/3)


def calculate_duration(period, rho, rprs, cosi):
    """
    Helper function to calculate transit duration predicted from a circular orbit
    
    Parameters
    ----------
    period : array-like
        orbital period [days]
    rho : float
        stellar density [solar density]
    rprs : float
        planet-to-star radius ratio for system
    cosi : array-like
        cosine(inclination)
        
    Returns
    -------
    transit_duration: array-like
        transit duration [days]
    """
    G = BIGG / RSUN**3 * MSUN * (24*3600)**2    # Newton's constant [R_sun^3 * M_sun^-1 * days^-2]
    
    term3  = ((3*period)/(G*rho*pi**2))**(1/3)
    term2a = (1+rprs)**2
    term2b = ((G*rho)/(3*pi))**(2/3)
    term2c = period**(4/3)*cosi**2
    term2  = (term2a - term2b*term2c)**(1/2)
    
    return term3*term2


def residuals_for_duration_fit(x0, x1, data_dur, data_err):
    """
    Helper function to return residuals for least squares fitting (op.leastsq)
    
    Parameters
    ----------
    x0 : array-like
        vector of parameters to vary in fit (cosi)
    x1 : array-like
        vector of parameters to hold constant (periods, rhostar, rprs)
    data_dur : array-like
        measured transit durations [days]
    data_err : array-like
        corresponding errors [days]
        
    Returns
    -------
    residuals : array-like
        error-scaled residuals on transit durations
    """
    cosi = x0
    period, rho, rprs = x1
    
    model_dur = calculate_duration(period, rho, rprs, cosi)
    
    return (data_dur - model_dur)/data_err


def calculate_flatness(data_dur, model_dur):
    """
    Helper function to calculate flatness
    
    Parameters
    ----------
    data_dur : array-like
        measured transit durations [days]
    model_dur : array-like
        model transit durations [days] from leastsq fit
        
    Returns
    -------
    flatness : array-like
        flatness measure
    """
    return np.std(data_dur-model_dur)/np.sqrt(np.mean(data_dur**2))



# Functions to compute system-level complexity measures
# Quantities definied in Gilbert & Fabrycky (2019)


def dynamical_mass(mp, Mstar):
    """
    Dynamical mass, mu
    
    Parameters
    ----------
    mp : array-like
        planet masses [M_earth]
    Mstar : float
        stellar mass [M_sun]
    """
    return np.sum(mp)/MSME/Mstar


def mass_partitioning(masses):
    """
    Mass partitioning, Q
    
    Parameters
    ----------
    masses : array-like
        planet masses [any units]
    """
    return LMC.D(masses/np.sum(masses))


def monotonicity(periods, masses):
    """
    Monotonicity, M
    
    Parameters
    ----------
    periods : array-like
        orbital periods [any units]
    masses : array-like
        planet masses corresponding to each given period [any units]
    """
    N = len(periods)
    rho = stats.spearmanr(periods, masses)[0]
    Q = mass_partitioning(masses)
    
    
    return rho*Q**(1/N)


def characteristic_spacing(periods, mp, Mstar, warn=True):
    """
    Characteristic spacing, S
    
    Parameters
    ----------
    periods : array-like
        orbital periods [days]
    mp : array-like
        planet masses corresponding to each given period [M_earth]
    Mstar : float
        Stellar mass [M_sun]
    warn : bool (optional)
        flag to control warnings (default=True)    
    """
    if len(periods) < 2:
        if warn:
            warnings.warn('Characteristic spacing is undefined for S < 2; returning NaN')
        return np.nan
    
    elif len(periods) >= 2:
        order = np.argsort(periods)
    
        periods = periods[order]
        mp = mp[order]
    
        a = P_to_a(periods, Mstar)

        radius_H = ((mp[1:]+mp[:-1])/(3*Mstar*MSME))**(1/3) * (a[1:]+a[:-1])/2
        delta_H  = (a[1:]-a[:-1])/radius_H

        return np.mean(delta_H)


def gap_complexity(periods, warn=True):
    """
    Gap complexity, C
    
    Parameters
    ----------
    periods : array-like
        planet periods
    warn : bool (optional)
        flag to control warnings (default=True)
    """
    if len(periods) < 3:
        if warn:
            warnings.warn('Complexity is undefined for N < 3; returning NaN')
        return np.nan
    
    elif len(periods) >= 3:
        order = np.argsort(periods)
  
        P = np.array(periods)[order]
        pp = np.log(P[1:]/P[:-1])/np.log(P.max()/P.min())
        
        return LMC.C(pp)
    

def flatness(periods, rhostar, rprs, dur, dur_err):
    """
    Flatness, f
    
    Parameters
    ----------
    periods : array-like
        orbital periods [days]
    rhostar : float
        stellar density [solar density]
    rprs : array-like
        planet-to-star radius ratios corresponding to given periods
    dur : array-like
        transit durations corresponding to given periods [hours]
    dur_err : array-like
        corresponding errors on transit durations [hours
    """
    cosi = np.array([0.0])
    transit_params = [periods, rhostar, rprs]
        
    cosi, success = op.leastsq(residuals_for_duration_fit, cosi, args=(transit_params, dur, dur_err))
        
    model_dur = calculate_duration(periods, rhostar, rprs, cosi)
    
    return calculate_flatness(dur, model_dur)
