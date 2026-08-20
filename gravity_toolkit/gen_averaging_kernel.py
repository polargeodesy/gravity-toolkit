#!/usr/bin/env python
"""
gen_averaging_kernel.py
Original IDL code gen_wclms_me.pro written by Sean Swenson
Adapted by Tyler Sutterley (06/2023)

Generates averaging kernel coefficients which minimize the total error

CALLING SEQUENCE:
    Wlms = gen_averaging_kernel(
        gclm, gslm, eclm, eslm, sigma,
        LMIN=0, LMAX=60, UNITS=0, RAD=300, LOVE=(hl,kl,ll)
    )

INPUTS:
    gclm: cosine spherical harmonics of exact averaging kernel
    gslm: sine spherical harmonics of exact averaging kernel
    eclm: measurement error in the cosine harmonics
    eslm: measurement error in the sine harmonics
    sigma: variance of the surface mass signal

OPTIONS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    RAD: Gaussian radius of the kernel (km)
    UNITS: units of input spherical harmonics
        0: fully-normalized
        1: mass coefficients (cmwe)
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

OUTPUTS:
    clm: cosine coefficients of the averaging kernel
    slm: sine coefficients of the averaging kernel

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    units.py: class for converting spherical harmonic data to specific units

REFERENCES:
    Swenson and Wahr, "Methods for inferring regional surface-mass anomalies
        from Gravity Recovery and Climate Experiment (GRACE) measurements of
        time-variable gravity," Journal of Geophysical Research: Solid Earth,
        107(B9), (2002). https://doi.org/10.1029/2001JB000576

UPDATE HISTORY:
    Updated 08/2026: separate the Gaussian kernel function
        rename Gaussian half-width to RAD to parallel other functions
    Updated 06/2023: added option for setting minimum value threshold
        use harmonics class for spherical harmonic operations
    Updated 04/2023: allow love numbers to be None for mass units case
    Updated 03/2023: improve typing for variables in docstrings
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2021: using units module for Earth parameters
    Updated 04/2020: reading load love numbers outside of this function
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Written 05/2013
"""

import numpy as np
from gravity_toolkit import gauss_kernel, harmonics, units


def gen_averaging_kernel(
    gclm,
    gslm,
    eclm,
    eslm,
    sigma,
    LMAX=60,
    MMAX=None,
    RAD=0,
    UNITS=0,
    LOVE=None,
    CUTOFF=1e-15,
):
    r"""
    Generates averaging kernel coefficients which minimize the
    total error following :cite:t:`Swenson:2002hs`

    Uses a normalized form of the Gaussian averaging function
    from :cite:p:`Jekeli:1981vj`

    Parameters
    ----------
    gclm: np.ndarray
        cosine spherical harmonics of exact averaging kernel
    gslm: np.ndarray
        sine spherical harmonics of exact averaging kernel
    eclm: np.ndarray
        measurement error in the cosine harmonics
    eslm: np.ndarray
        measurement error in the sine harmonics
    sigma: float
        variance of the surface mass signal
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    RAD: float, default 0
        Gaussian radius of the kernel (km)
    UNITS: int, default 0
        Input data units

            - ``0``: fully-normalized
            - ``1``: mass coefficients (cm w.e., g/cm\ :sup:`2`)
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    CUTOFF: float, default 1e-15
        minimum value for tail of Gaussian kernel

    Returns
    -------
    clm: np.ndarray
        cosine coefficients of the averaging kernel
    slm: np.ndarray
        sine coefficients of the averaging kernel
    """
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # Earth Parameters
    factors = units(lmax=LMAX)
    # extract arrays of kl, hl, and ll Love Numbers
    if UNITS == 0:
        # Input coefficients are fully-normalized
        dfactor = factors.harmonic(*LOVE).cmwe
    elif UNITS == 1:
        # Inputs coefficients are mass (cmwe)
        dfactor = np.ones((LMAX + 1))

    # calculate legendre coefficients of a Gaussian correlation function
    gl = gauss_kernel(RAD, LMAX, CUTOFF=CUTOFF)

    # copy of the area under the kernel
    area = np.copy(gclm[0, 0])

    # Convert sigma to correlation function amplitude
    variance = np.zeros((LMAX + 1))
    for l in range(0, LMAX + 1):  # equivalent to 0:LMAX
        mm = np.min([MMAX, l])  # find min of MMAX and l
        m = np.arange(0, mm + 1)  # create m array 0:l or 0:MMAX
        variance[l] = (gl[l] / 2.0) * np.sum(gclm[l, m] ** 2 + gslm[l, m] ** 2)

    # divide by the square of the area under the kernel
    root = np.sqrt(np.sum(variance) / np.power(area, 2))
    # signal variance
    sigma_0 = sigma / root

    # Compute averaging kernel coefficients
    Ylms = harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX + 1, MMAX + 1))
    Ylms.slm = np.zeros((LMAX + 1, MMAX + 1))
    # for each spherical harmonic degree
    for l in range(0, LMAX + 1):  # equivalent to 0:lmax
        # inverse of smoothed signal variance in output units
        coeff = (dfactor[l] ** 2) / (gl[l] * sigma_0**2)
        # for each valid spherical harmonic order
        mm = np.min([MMAX, l])
        m = np.arange(0, mm + 1)
        # compute averaging kernel coefficients
        temp = 1.0 + 2.0 * coeff * eclm[l, m] ** 2
        Ylms.clm[l, m] = gclm[l, m] / temp
        temp = 1.0 + 2.0 * coeff * eslm[l, m] ** 2
        Ylms.slm[l, m] = gslm[l, m] / temp

    # return kernels divided by the area under the kernel
    return Ylms.scale(1.0 / area)
