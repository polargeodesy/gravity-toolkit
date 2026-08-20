=================================
``gravity_toolkit.gauss_weights``
=================================

- Computes the Gaussian weights as a function of degree
- Normalized form of the Gaussian averaging function from :cite:t:`Jekeli:1981vj`

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gauss_weights import gauss_weights
    wl = 2.0*np.pi*gauss_weights(hw,LMAX)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/gauss_weights.py

.. autofunction:: gravity_toolkit.gauss_weights

.. autofunction:: gravity_toolkit.gauss_kernel
