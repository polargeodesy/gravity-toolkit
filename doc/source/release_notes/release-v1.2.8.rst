.. _release-v1.2.8:

==================
`Release v1.2.8`__
==================

* ``refactor``: adjust format of logfiles `(#187) <https://github.com/polargeodesy/gravity-toolkit/pull/187>`_
* ``refactor``: reorganize for modern ``setuptools`` builds `(#187) <https://github.com/polargeodesy/gravity-toolkit/pull/187>`_
* ``fix``: simplify Clenshaw summation to reduce memory usage `(#187) <https://github.com/polargeodesy/gravity-toolkit/pull/187>`_
* ``docs``: fix some broken links `(#189) <https://github.com/polargeodesy/gravity-toolkit/pull/189>`_
* ``feat``: save start and end date of files `(#190) <https://github.com/polargeodesy/gravity-toolkit/pull/190>`_
* ``refactor``: output rasters using structured netCDF4 function from ``geoid-toolkit`` `(#190) <https://github.com/polargeodesy/gravity-toolkit/pull/190>`_
* ``fix``: include additional attributes to output files for CF compliance `(#190) <https://github.com/polargeodesy/gravity-toolkit/pull/190>`_
* ``refactor``: spherical harmonic errors using euler's formula and ``np.einsum`` `(#190) <https://github.com/polargeodesy/gravity-toolkit/pull/190>`_
* ``fix``: add ``astype`` to ``calendar_days`` function `(#190) <https://github.com/polargeodesy/gravity-toolkit/pull/190>`_
* ``refactor``: add ocean and land geocenter functions to class `(#192) <https://github.com/polargeodesy/gravity-toolkit/pull/192>`_
* ``refactor``: change assertions to value errors `(#192) <https://github.com/polargeodesy/gravity-toolkit/pull/192>`_
* ``test``: add geocenter tests `(#192) <https://github.com/polargeodesy/gravity-toolkit/pull/192>`_
* ``refactor``: separate the Gaussian kernel function `(#194) <https://github.com/polargeodesy/gravity-toolkit/pull/194>`_
* ``test``: add a gaussian weight test `(#194) <https://github.com/polargeodesy/gravity-toolkit/pull/194>`_
* ``chore``: transfer ownership to ``polargeodesy`` `(#195) <https://github.com/polargeodesy/gravity-toolkit/pull/195>`_
* ``fix``: allocate for and then fill output spherical harmonics `(#196) <https://github.com/polargeodesy/gravity-toolkit/pull/196>`_
* ``fix``: check dimensions of input element x `(#196) <https://github.com/polargeodesy/gravity-toolkit/pull/196>`_
* ``fix``: ``plt.get_cmap`` instead of ``cm.get_cmap`` `(#196) <https://github.com/polargeodesy/gravity-toolkit/pull/196>`_
* ``refactor``: format ``pyproject.toml`` with ``taplo`` `(#196) <https://github.com/polargeodesy/gravity-toolkit/pull/196>`_
* ``fix``: root attributes for spatial output `(#197) <https://github.com/polargeodesy/gravity-toolkit/pull/197>`_
* ``fix``: convert back to just using ``np.sum`` for point loads `(#199) <https://github.com/polargeodesy/gravity-toolkit/pull/199>`_
* ``feat``: add options to use ridge regression with tunable lambda `(#201) <https://github.com/polargeodesy/gravity-toolkit/pull/201>`_

.. __: https://github.com/polargeodesy/gravity-toolkit/releases/tag/1.2.8
