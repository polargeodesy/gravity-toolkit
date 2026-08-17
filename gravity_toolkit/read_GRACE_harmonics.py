#!/usr/bin/env python
"""
read_GRACE_harmonics.py
Written by Tyler Sutterley (08/2026)
Contributions by Hugo Lecomte

Reads GRACE files and extracts spherical harmonic data and drift rates (RL04)
Adds drift rates to clm and slm for release 4 harmonics
Correct GSM data for drift in pole tide following Wahr et al. (2015)
Parses date of GRACE/GRACE-FO data from filename

INPUTS:
    input_file: GRACE/GRACE-FO Level-2 spherical harmonic data file
    LMAX: Maximum degree of spherical harmonics (degree of truncation)

OPTIONS:
    MMAX: Maximum order of spherical harmonics (order of truncation)
        default is the maximum spherical harmonic degree
    POLE_TIDE: correct GSM data for pole tide drift following Wahr et al. (2015)

OUTPUTS:
    time: mid-month date in year-decimal
    start: start date of range as Julian day
    end: end date of range as Julian day
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic uncalibrated standard deviations
    eslm: sine spherical harmonic uncalibrated standard deviations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    PyYAML: YAML parser and emitter for Python
        https://github.com/yaml/pyyaml

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 08/2026: use python datetime to calculate start and end dates
    Updated 11/2024: check if the GRACE/GRACE-FO files are gfc format
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: added regex formatting for CNES GRGS harmonics
        improve typing for variables in docstrings
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: make keyword arguments part of kwargs dictionary
    Updated 05/2022: updated comments
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
        check if GRACE/GRACE-FO data file is present in file-system
    Updated 09/2021: added COST-G combined solutions from the GFZ ICGEM
        output spherical harmonic degree and order in dict
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 12/2020: using utilities from time module
    Updated 10/2020: Change parse function to work with GRGS data
    Updated 08/2020: flake8 compatible regular expression strings
        input file can be "diskless" bytesIO object
    Updated 07/2020: added function docstrings
    Updated 08/2019: specify yaml loader (PyYAML yaml.load(input) Deprecation)
    Updated 07/2019: replace colons in yaml header if within quotations
    Updated 11/2018: decode gzip read with ISO-8859-1 for python3 compatibility
    Updated 05/2018: updates to file name structure with release 6 and GRACE-FO
        output file headers and parse new YAML headers for RL06 and GRACE-FO
    Written 10/2017 for public release
"""

import re
import io
import gzip
import yaml
import pathlib
import numpy as np
import gravity_toolkit.time
from datetime import datetime, timedelta


# PURPOSE: read Level-2 GRACE and GRACE-FO spherical harmonic files
def read_GRACE_harmonics(input_file, LMAX, **kwargs):
    """
    Extracts spherical harmonic coefficients from GRACE/GRACE-FO files

    Parameters
    ----------
    input_file: str
        GRACE/GRACE-FO Level-2 spherical harmonic data file
    LMAX: int
        Maximum degree of spherical harmonics (degree of truncation)
    MMAX: int or NoneType, default None
        Maximum order of spherical harmonics
    POLE_TIDE: bool, default False
        Correct for pole tide drift following :cite:t:`Wahr:2015dg`

    Returns
    -------
    time: float
        mid-month date in year-decimal
    start: float
        start date of range as Julian day
    end: float
        end date of range as Julian day
    l: np.ndarray
        spherical harmonic degree to LMAX
    m: np.ndarray
        spherical harmonic order to MMAX
    clm: np.ndarray
        cosine spherical harmonics coefficients
    slm: np.ndarray
        sine spherical harmonics coefficients
    eclm: np.ndarray
        cosine spherical harmonic uncalibrated standard deviations
    eslm: np.ndarray
        sine spherical harmonic uncalibrated standard deviations
    header: str
        Header text from the GRACE/GRACE-FO file
    """
    # set default keyword arguments
    kwargs.setdefault('MMAX', None)
    kwargs.setdefault('POLE_TIDE', False)

    # parse filename
    PFX, SY, SD, EY, ED, N, PRC, F1, DRL, F2, SFX = parse_file(input_file)
    # check if file is compressed
    compressed = SFX == '.gz'
    # extract file contents
    file_contents = extract_file(input_file, compressed)

    # JPL mascon solutions in spherical harmonic form
    if PRC in ('JPLMSC',):
        DSET = 'GSM'
        DREL = np.int64(DRL)
        FLAG = r'GRCOF2'
    # Kusche et al. (2009) DDK filtered solutions
    # https://doi.org/10.1007/s00190-009-0308-3
    elif PFX.startswith('kfilter_DDK'):
        DSET = 'GSM'
        DREL = np.int64(DRL)
        FLAG = r'gfc'
    # COST-G unfiltered combination solutions
    # https://doi.org/10.5880/ICGEM.COST-G.001
    # GFC solutions from the GFZ ICGEM
    # https://icgem.gfz-potsdam.de/sl/temporal
    elif PRC in ('COSTG',) or SFX in ('.gfc',):
        (DSET,) = re.findall(r'(GSM|GAA|GAB|GAC|GAD)', PFX)
        DREL = np.int64(DRL)
        FLAG = r'gfc'
    # Standard GRACE/GRACE-FO Level-2 solutions
    else:
        (DSET,) = re.findall(r'(GSM|GAA|GAB|GAC|GAD)', PFX)
        DREL = np.int64(DRL)
        FLAG = r'GRCOF2'

    # output python dictionary with GRACE/GRACE-FO data and metadata
    # spherical harmonic model (SHM) data
    SHM = {}

    # extract GRACE/GRACE-FO date information from input file name
    start_date = datetime(int(SY), 1, 1) + timedelta(days=int(SD) - 1)
    start_struct = start_date.timetuple()
    end_date = datetime(int(EY), 1, 1) + timedelta(days=int(ED) - 1)
    end_struct = end_date.timetuple()
    # start and end day of the year
    start_yr = start_struct.tm_year
    start_day = start_struct.tm_yday
    end_yr = end_struct.tm_year
    end_day = end_struct.tm_yday
    # calculate mid-month date taking into account if measurements are
    # on different years
    dpy = gravity_toolkit.time.calendar_days(start_yr).sum()

    # For data that crosses years (end_yr - start_yr should be at most 1)
    end_cyclic = (end_yr - start_yr) * dpy + end_day
    # Calculate mid-month value
    mid_day = np.mean([start_day, end_cyclic])
    # Calculating the mid-month date in decimal form
    SHM['time'] = start_yr + mid_day / dpy

    # Calculating the Julian dates of the start and end date
    MJD1 = gravity_toolkit.time.convert_calendar_dates(
        start_yr,
        start_struct.tm_mon,
        start_struct.tm_mday,
        epoch=(1858, 11, 17, 0, 0, 0),
    )
    MJD2 = gravity_toolkit.time.convert_calendar_dates(
        end_yr,
        end_struct.tm_mon,
        end_struct.tm_mday,
        epoch=(1858, 11, 17, 0, 0, 0),
    )
    SHM['start'] = 2400000.5 + MJD1
    SHM['end'] = 2400000.5 + MJD2

    # set maximum spherical harmonic order
    MMAX = kwargs.get('MMAX', None)
    # only replace if None (allow MMAX to be zero, which is typically falsy)
    if MMAX is None:
        MMAX = np.copy(LMAX)
    # output dimensions
    SHM['l'] = np.arange(LMAX + 1)
    SHM['m'] = np.arange(MMAX + 1)
    # Spherical harmonic coefficient matrices to be filled from data file
    SHM['clm'] = np.zeros((LMAX + 1, MMAX + 1))
    SHM['slm'] = np.zeros((LMAX + 1, MMAX + 1))
    # spherical harmonic uncalibrated standard deviations
    SHM['eclm'] = np.zeros((LMAX + 1, MMAX + 1))
    SHM['eslm'] = np.zeros((LMAX + 1, MMAX + 1))
    if (DREL == 4) and (DSET == 'GSM'):
        # clm and slm drift rates for RL04
        drift_c = np.zeros((LMAX + 1, MMAX + 1))
        drift_s = np.zeros((LMAX + 1, MMAX + 1))
    # set default degree 0 harmonics for intercomparability between centers
    SHM['clm'][0, 0] = 1.0

    # extract GRACE and GRACE-FO file headers
    # replace colons in header if within quotations
    head = [
        re.sub(r'\"(.*?)\:\s(.*?)\"', r'"\1, \2"', l)
        for l in file_contents
        if not re.match(rf'{FLAG}|GRDOTA', l)
    ]
    if SFX in ('.gfc',):
        # extract parameters from header
        header_parameters = [
            'modelname',
            'earth_gravity_constant',
            'radius',
            'max_degree',
            'errors',
            'norm',
            'tide_system',
        ]
        header_regex = re.compile(r'(' + r'|'.join(header_parameters) + r')')
        header = [l.split(maxsplit=1) for l in head if header_regex.match(l)]
        SHM['header'] = {i[0]: i[1] for i in header}
    elif ((N == 'GRAC') and (DREL >= 6)) or (N == 'GRFO'):
        # parse the YAML header for RL06 or GRACE-FO (specifying yaml loader)
        SHM.update(yaml.load('\n'.join(head), Loader=yaml.BaseLoader))
    else:
        # save lines of the GRACE file header removing empty lines
        SHM['header'] = [l.rstrip() for l in head if l]

    # for each line in the GRACE/GRACE-FO file
    for line in file_contents:
        # find if line starts with data marker flag (e.g. GRCOF2)
        if bool(re.match(FLAG, line)):
            # split the line into individual components
            line_contents = line.split()
            # degree and order for the line
            l1 = np.int64(line_contents[1])
            m1 = np.int64(line_contents[2])
            # if degree and order are below the truncation limits
            if (l1 <= LMAX) and (m1 <= MMAX):
                SHM['clm'][l1, m1] = np.float64(line_contents[3])
                SHM['slm'][l1, m1] = np.float64(line_contents[4])
                SHM['eclm'][l1, m1] = np.float64(line_contents[5])
                SHM['eslm'][l1, m1] = np.float64(line_contents[6])
        # find if line starts with drift rate flag
        elif bool(re.match(r'GRDOTA', line)):
            # split the line into individual components
            line_contents = line.split()
            l1 = np.int64(line_contents[1])
            m1 = np.int64(line_contents[2])
            # Reading Drift rates for low degree harmonics
            drift_c[l1, m1] = np.float64(line_contents[3])
            drift_s[l1, m1] = np.float64(line_contents[4])

    # Adding drift rates to clm and slm for RL04
    # if drift rates exist at any time, will add to harmonics
    # Will convert the secular rates into a stokes contribution
    # Currently removes 2003.3 to get the temporal average close to 0.
    if (DREL == 4) and (DSET == 'GSM'):
        # time since 2003.3
        dt = SHM['time'] - 2003.3
        SHM['clm'][:, :] += dt * drift_c[:, :]
        SHM['slm'][:, :] += dt * drift_s[:, :]

    # Correct Pole Tide following Wahr et al. (2015) 10.1002/2015JB011986
    if kwargs['POLE_TIDE'] and (DSET == 'GSM'):
        # time since 2000.0
        dt = SHM['time'] - 2000.0
        # CSR and JPL Pole Tide Correction
        if PRC in ('UTCSR', 'JPLEM', 'JPLMSC'):
            # values for IERS mean pole [2010]
            if SHM['time'] < 2010.0:
                a = np.array([0.055974, 1.8243e-3, 1.8413e-4, 7.024e-6])
                b = np.array([-0.346346, -1.7896e-3, 1.0729e-4, 0.908e-6])
            elif SHM['time'] >= 2010.0:
                a = np.array([0.023513, 7.6141e-3, 0.0, 0.0])
                b = np.array([-0.358891, 0.6287e-3, 0.0, 0.0])
            # calculate m1 and m2 values
            m1 = np.copy(a[0])
            m2 = np.copy(b[0])
            for x in range(1, 4):
                m1 += a[x] * dt**x
                m2 += b[x] * dt**x
            # pole tide values for CSR and JPL
            # CSR and JPL both remove the IERS mean pole from m1 and m2
            # before computing their harmonic solutions
            C21_PT = -1.551e-9 * (m1 - 0.62e-3 * dt) - 0.012e-9 * (
                m2 + 3.48e-3 * dt
            )
            S21_PT = 0.021e-9 * (m1 - 0.62e-3 * dt) - 1.505e-9 * (
                m2 + 3.48e-3 * dt
            )
            # correct GRACE/GRACE-FO spherical harmonics for pole tide
            SHM['clm'][2, 1] -= C21_PT
            SHM['slm'][2, 1] -= S21_PT
        # GFZ Pole Tide Correction
        elif PRC in ('EIGEN', 'GFZOP'):
            # pole tide values for GFZ
            # GFZ removes only a constant pole position
            C21_PT = -1.551e-9 * (-0.62e-3 * dt) - 0.012e-9 * (3.48e-3 * dt)
            S21_PT = 0.021e-9 * (-0.62e-3 * dt) - 1.505e-9 * (3.48e-3 * dt)
            # correct GRACE/GRACE-FO spherical harmonics for pole tide
            SHM['clm'][2, 1] -= C21_PT
            SHM['slm'][2, 1] -= S21_PT

    # return the header data, GRACE/GRACE-FO data
    # GRACE/GRACE-FO date (mid-month in decimal)
    # and the start and end days as Julian dates
    return SHM


# PURPOSE: extract parameters from filename
def parse_file(input_file):
    """
    Extract parameters from filename

    Parameters
    ----------
    input_file: str
        GRACE/GRACE-FO Level-2 spherical harmonic data file
    """
    # compile numerical expression operator for parameters from files
    # UTCSR: The University of Texas at Austin Center for Space Research
    # EIGEN: GFZ German Research Center for Geosciences (RL01-RL05)
    # GFZOP: GFZ German Research Center for Geosciences (RL06+GRACE-FO)
    # JPLEM: NASA Jet Propulsion Laboratory (harmonic solutions)
    # JPLMSC: NASA Jet Propulsion Laboratory (mascon solutions)
    # GRGS: French Centre National D'Etudes Spatiales (CNES)
    # COSTG: International Combined Time-variable Gravity Fields
    # GRGS: CNES Groupe de Recherche de Geodesie Spatiale
    centers = r'UTCSR|EIGEN|GFZOP|JPLEM|JPLMSC|GRGS|COSTG|GRGS'
    suffixes = r'\.gz|\.gfc|\.txt'
    regex_pattern = (
        r'(.*?)-2_(\d{4})(\d{3})-(\d{4})(\d{3})_'
        rf'(.*?)_({centers})_(.*?)_(\d+)(.*?)({suffixes})?$'
    )
    rx = re.compile(regex_pattern, re.VERBOSE)
    # extract parameters from input filename
    if isinstance(input_file, io.IOBase):
        return rx.findall(input_file.filename).pop()
    else:
        return rx.findall(pathlib.Path(input_file).name).pop()


# PURPOSE: read input file and extract contents
def extract_file(input_file, compressed):
    """
    Read input file and extract contents

    Parameters
    ----------
    input_file: str
        GRACE/GRACE-FO Level-2 spherical harmonic data file
    compressed: bool
        denotes if the file is compressed
    """
    # tilde expansion of input file if not byteIO object
    if not isinstance(input_file, io.IOBase):
        input_file = pathlib.Path(input_file).expanduser().absolute()
        # check that data file is present in file system
        if not input_file.exists():
            raise FileNotFoundError(f'{str(input_file)} not found')
    # check if file is uncompressed byteIO object
    if isinstance(input_file, io.IOBase) and not compressed:
        # extract spherical harmonic coefficients
        return input_file.read().decode('ISO-8859-1').splitlines()
    else:
        # check if file is compressed (read with gzip if gz)
        file_opener = gzip.open if compressed else open
        # opening data file to extract spherical harmonic coefficients
        with file_opener(input_file, 'rb') as f:
            return f.read().decode('ISO-8859-1').splitlines()
