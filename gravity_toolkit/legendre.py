#!/usr/bin/env python
"""
legendre.py
Written by Tyler Sutterley (09/2026)
Computes associated Legendre functions of degree l evaluated for elements x
l must be a scalar integer and x must contain real values ranging -1 <= x <= 1
Parallels the MATLAB legendre function

Based on Fortran program by Robert L. Parker, Scripps Institution of
Oceanography, Institute for Geophysics and Planetary Physics, UCSD. 1993

INPUTS:
    l: degree of Legendre polynomials
    x: elements ranging from -1 to 1
        typically cos(theta), where theta is the colatitude in radians

OUTPUT:
    Pl: legendre polynomials of degree l for orders 0 to l

OPTIONS:
    NORMALIZE: output Fully Normalized Associated Legendre Functions

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

REFERENCES:
    M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
        Dover Publications, 1965, Ch. 8.
    J. A. Jacobs, "Geomagnetism", Academic Press, 1987, Ch.4.

UPDATE HISTORY:
    Updated 09/2026: check dimensions of input element x
    Updated 08/2026: fixes for typing error with numpy updates
    Updated 03/2023: improve typing for variables in docstrings
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: modify normalization to prevent high degree overflows
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: modify case with underflow
    Updated 09/2020: verify dimensions of x variable
    Updated 07/2020: added function docstrings
    Updated 05/2020: added normalization option for output polynomials
    Updated 03/2019: calculate twocot separately to avoid divide warning
    Written 08/2016
"""

import numpy as np


def legendre(l, x, NORMALIZE=False):
    r"""
    Computes associated Legendre functions for a particular degree
    following :cite:t:`Abramowitz:1965vw,Jacobs:1987vv`

    Parameters
    ----------
    l: int
        degree of Legendre polynomials
    x: np.ndarray
        elements ranging from -1 to 1

        Typically :math:`\cos(\theta)`, where :math:`\theta`
        is the colatitude in radians
    NORMALIZE: bool, default False
        Fully-normalize the Legendre Functions

    Returns
    -------
    Pl: np.ndarray
        legendre polynomials of degree ``l``
    """
    # verify integer
    l = np.int64(l)
    # check dimensions of input elements
    singular_values = np.ndim(x) == 0
    # verify dimensions
    x = np.atleast_1d(x).flatten()
    # size of the x array
    nx = len(x)
    # tolerances for underflow
    tol = np.sqrt(np.finfo(np.float64).tiny)
    tstart = np.finfo(np.float64).eps

    # for the l = 0 case
    if l == 0:
        Pl = np.ones((1, nx), dtype=np.float64)
        return Pl

    # for all other degrees greater than 0
    rootl = np.sqrt(np.arange(0, 2 * l + 1))  # +1 to include 2*l
    # u is sine of colatitude (cosine of latitude) so that 0 <= u <= 1
    u = np.sqrt(1.0 - x**2)
    P = np.zeros((l + 3, nx), dtype=np.float64)
    # sine of -colatitude to power l
    upow = np.power(-u, l)

    # calculate legendre polynomials for underflow cases
    if np.any((u > 0) & (np.abs(upow) <= tol)):
        # find indices where the power terms are small
        ind = np.flatnonzero((u > 0) & (np.abs(upow) <= tol))
        # approximate solution of x*ln(x) = Pl
        v = 9.2 - np.log(tol) / (l * u[ind])
        w = 1.0 / np.log(v)
        m1 = 1 + l * u[ind] * v * w * (1.0058 + w * (3.819 - w * 12.173))
        m1 = np.where(l < np.floor(m1), l, np.floor(m1)).astype(np.int64)
        # column-by-column recursion
        for k, mm1 in enumerate(m1):
            col = ind[k]
            # calculate two*cotangent for underflow case
            twocot = -2.0 * x[col] / u[col]
            P[mm1 - 1 : l + 1, col] = 0.0
            # start recursion with proper sign
            P[mm1 - 1, col] = np.sign(np.fmod(mm1, 2) - 0.5) * tstart
            if x[col] < 0:
                P[mm1 - 1, col] = np.sign(np.fmod(l + 1, 2) - 0.5) * tstart
            # backwards recursion from m1 to m = 0
            # accumulate the normalization factor
            sumsq = tol.copy()
            for m in range(mm1 - 2, -1, -1):
                P[m, col] = (
                    (m + 1) * twocot * P[m + 1, col]
                    - rootl[l + m + 2] * rootl[l - m - 1] * P[m + 2, col]
                ) / (rootl[l + m + 1] * rootl[l - m])
                sumsq += P[m, col] ** 2
            # calculate scale
            scale = 1.0 / np.sqrt(2.0 * sumsq - P[0, col] ** 2)
            P[0 : mm1 + 1, col] = scale * P[0 : mm1 + 1, col]

    # calculate legendre polynomials for most (normal) cases
    # (no underflow and not polar)
    if np.any((x != 1) & (np.abs(upow) >= tol)):
        # find indices where the power terms are above tolerance
        ind = np.flatnonzero((x != 1) & (np.abs(upow) >= tol))
        # calculate two*cotangent for normal case
        twocot = -2.0 * x[ind] / u[ind]
        # produce normalization constant for the m = l function
        d = np.arange(2, 2 * l + 2, 2)
        c = np.prod(1.0 - 1.0 / d)
        # calculate for degree l and use backwards recursion for other degrees
        P[l, ind] = np.sqrt(c) * upow[ind]
        P[l - 1, ind] = P[l, ind] * twocot * l / rootl[-1]
        # recur downwards to m = 0
        for m in range(l - 2, -1, -1):
            P[m, ind] = (
                P[m + 1, ind] * twocot * (m + 1)
                - P[m + 2, ind] * rootl[l + m + 2] * rootl[l - m - 1]
            ) / (rootl[l + m + 1] * rootl[l - m])

    # polar case (x == +/-1)
    if np.any(u == 0):
        # find indices for values at the poles
        u0 = np.flatnonzero(u == 0)
        P[0, u0] = np.power(x[u0], l)

    # calculate Pl from P
    # and truncate to degree l
    Pl = np.copy(P[0 : l + 1, :])

    # calculate Fully Normalized Associated Legendre functions
    if NORMALIZE:
        # allocate for normalization array
        norm = np.zeros((l + 1))
        # normalization for degree 0
        norm[0] = np.sqrt(2.0 * l + 1)
        # normalization for all other degrees
        m = np.arange(1, l + 1)
        # apply the Condon-Shortley phase
        cs = np.power(-1.0, m)
        norm[1:] = cs * np.sqrt(2.0 * (2.0 * l + 1.0))
        # apply normalization to each row of Pl
        # reshape the normalization array to be compatible with Pl
        Pl *= np.kron(np.ones((1, nx)), norm[:, np.newaxis])
    else:
        # Calculate the unnormalized Legendre functions by multiplying each row
        # by: sqrt((l+m)!/(l-m)!) == sqrt(prod(n-m+1:n+m))
        # following Abramowitz and Stegun
        for m in range(1, l):
            Pl[m, :] *= np.prod(rootl[l - m + 1 : l + m + 1])
        # sectoral case (l = m) should be done separately to handle 0!
        Pl[l, :] *= np.prod(rootl[1:])

    # return the legendre polynomials
    # flattened to singular values if necessary
    if singular_values:
        return Pl[:, 0]
    else:
        return Pl
