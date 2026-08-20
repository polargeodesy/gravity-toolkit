==============================
``dealiasing_monthly_mean.py``
==============================

- Reads GRACE/GRACE-FO level-1b dealiasing data files for a specific product and outputs monthly the mean for a specific GRACE/GRACE-FO processing center and data release

    * ``'GAA'``: atmospheric loading from ECMWF
    * ``'GAB'``: oceanic loading from OMCT/MPIOM
    * ``'GAC'``: global atmospheric and oceanic loading
    * ``'GAD'``: ocean bottom pressure from OMCT/MPIOM
- Creates monthly files of oblateness variations at 3 or 6-hour intervals

`Source code`__

.. __: https://github.com/polargeodesy/gravity-toolkit/blob/main/gravity_toolkit/scripts/dealiasing_monthly_mean.py

Calling Sequence
################

.. argparse::
    :module: gravity_toolkit.scripts.dealiasing_monthly_mean
    :func: arguments
    :prog: dealiasing_monthly_mean.py
    :nodescription:
    :nodefault:
