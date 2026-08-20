==========================
``gfz_isdc_grace_sync.py``
==========================

- Syncs GRACE/GRACE-FO and auxiliary data from the `GFZ Information System and Data Center (ISDC) <http://isdc.gfz-potsdam.de/grace-isdc/>`_
- Syncs CSR/GFZ/JPL Level-2 spherical harmonic files
- Gets the latest technical note (TN) files
- Gets the monthly GRACE/GRACE-FO newsletters
- Creates an index file for each data product

`Source code`__

.. __: https://github.com/polargeodesy/gravity-toolkit/blob/main/gravity_toolkit/datasets/gfz_isdc_grace_sync.py

Calling Sequence
################

.. argparse::
    :module: gravity_toolkit.datasets.gfz_isdc_grace_sync
    :func: arguments
    :prog: gfz_isdc_grace_sync.py
    :nodescription:
    :nodefault:
