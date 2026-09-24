#!/usr/bin/env python
"""
test_gaussian.py (08/2026)
Tests Gaussian weights versus outputs from recurs.f
"""

import inspect
import pathlib
import numpy as np
import gravity_toolkit as gravtk

# path to test files
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent


# PURPOSE: test gaussian weights functions
def test_gaussian_weights():
    # read input Gaussian weights file
    filename = filepath.joinpath(f'out.recurs.csv.gz')
    table = np.genfromtxt(filename, delimiter=',', names=True)
    # maximum degree and order of table
    lmax = table['l'].max().astype('i')
    # table has four columns of values for alphas 450, 500, 625, 900
    for i, alpha in enumerate([450, 500, 625, 900]):
        # data column in the csv
        column = f'b_alpha_{alpha:d}'
        # convert from alpha to Gaussian half-width
        hw = 6371.0 * np.arccos(1.0 - np.log(2.0) / alpha)
        # calculate Gaussian weights
        # normalize to compare with outputs from recurs
        gw = 2.0 * np.pi * gravtk.gauss_weights(hw, lmax, CUTOFF=1e-10)
        # check validity
        assert np.allclose(gw, table[column], atol=5e-7)
