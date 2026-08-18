#!/usr/bin/env python
"""
test_geocenter.py (08/2026)
"""

import pytest
import numpy as np
import gravity_toolkit as gravtk


# PURPOSE: model the seasonal component of an initial degree 1 model
# using preliminary estimates of annual and semi-annual variations from LWM
# as calculated in Chen et al. (1999), doi:10.1029/1998JB900019
# NOTE: this is to get an accurate assessment of the land water mass for the
# eustatic component (not for the ocean component from GRACE)
def test_land_seasonal():
    # create a range of test dates
    grace_date = np.arange(2002.25, 2020.25, 1.0 / 12.0)
    # Annual amplitudes of (Soil Moisture + Snow) geocenter components (mm)
    AAx = 1.28
    AAy = 0.52
    AAz = 3.30
    # Annual phase of (Soil Moisture + Snow) geocenter components (degrees)
    APx = 44.0
    APy = 182.0
    APz = 43.0
    # Semi-Annual amplitudes of (Soil Moisture + Snow) geocenter components
    SAAx = 0.15
    SAAy = 0.56
    SAAz = 0.50
    # Semi-Annual phase of (Soil Moisture + Snow) geocenter components
    SAPx = 331.0
    SAPy = 312.0
    SAPz = 75.0
    # calculate each geocenter component from the amplitude and phase
    # converting the phase from degrees to radians
    X = AAx * np.sin(
        2.0 * np.pi * grace_date + np.radians(APx)
    ) + SAAx * np.sin(4.0 * np.pi * grace_date + np.radians(SAPx))
    Y = AAy * np.sin(
        2.0 * np.pi * grace_date + np.radians(APy)
    ) + SAAy * np.sin(4.0 * np.pi * grace_date + np.radians(SAPy))
    Z = AAz * np.sin(
        2.0 * np.pi * grace_date + np.radians(APz)
    ) + SAAz * np.sin(4.0 * np.pi * grace_date + np.radians(SAPz))
    valid = gravtk.geocenter(X=X - X.mean(), Y=Y - Y.mean(), Z=Z - Z.mean())
    valid.from_cartesian()
    # calculate using direct function
    DEG1 = gravtk.geocenter.land_seasonal(grace_date)
    # compare geocenter and degree one components
    for key in ['X', 'Y', 'Z', 'C10', 'C11', 'S11']:
        assert np.allclose(valid[key], DEG1[key])


def test_ocean_seasonal():
    # create a range of test dates
    grace_date = np.arange(2002.25, 2020.25, 1.0 / 12.0)
    # Annual amplitudes of ocean (TOPEX) geocenter components (mm)
    AAx = 0.96
    AAy = 0.97
    AAz = 0.49
    # Annual phase of ocean (TOPEX) geocenter components (degrees)
    APx = 73.0
    APy = 52.0
    APz = 3.0
    # Semi-Annual amplitudes of ocean (TOPEX) geocenter components
    SAAx = 0.86
    SAAy = 0.73
    SAAz = 0.25
    # Semi-Annual phase of ocean (TOPEX) geocenter components
    SAPx = 187.0
    SAPy = 173.0
    SAPz = 232.0
    # calculate each geocenter component from the amplitude and phase
    # converting the phase from degrees to radians
    X = AAx * np.sin(
        2.0 * np.pi * grace_date + np.radians(APx)
    ) + SAAx * np.sin(4.0 * np.pi * grace_date + np.radians(SAPx))
    Y = AAy * np.sin(
        2.0 * np.pi * grace_date + np.radians(APy)
    ) + SAAy * np.sin(4.0 * np.pi * grace_date + np.radians(SAPy))
    Z = AAz * np.sin(
        2.0 * np.pi * grace_date + np.radians(APz)
    ) + SAAz * np.sin(4.0 * np.pi * grace_date + np.radians(SAPz))
    valid = gravtk.geocenter(X=X - X.mean(), Y=Y - Y.mean(), Z=Z - Z.mean())
    valid.from_cartesian()
    # calculate using direct function
    DEG1 = gravtk.geocenter.ocean_seasonal(grace_date)
    # compare geocenter and degree one components
    for key in ['X', 'Y', 'Z', 'C10', 'C11', 'S11']:
        assert np.allclose(valid[key], DEG1[key])


# PURPOSE: test the from_harmonics class method
def test_from_harmonics():
    rng = np.random.default_rng()
    Ylms = gravtk.harmonics(lmax=1).zeros()
    Ylms.clm[1, 0] = rng.random()
    Ylms.clm[1, 1] = rng.random()
    Ylms.slm[1, 1] = rng.random()
    geocenter = gravtk.geocenter.from_harmonics(Ylms)
    assert geocenter.C10 == Ylms.clm[1, 0]
    assert geocenter.C11 == Ylms.clm[1, 1]
    assert geocenter.S11 == Ylms.slm[1, 1]
