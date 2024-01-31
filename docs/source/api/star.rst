.. module:: standard_classification

.. _api:

The star class
===============

This class is the main data structure of this package. it represents an
initially unclassified YSO to which the standard classification is
applied to. To create a new instance of the star class a astropy 
:code:`Table.Row` object is used. It holds the photometric data for the
given source. The data should consist of spectral flux densities given
in Jansky (:math:`[F_\nu] = \mathrm{Jy}`) and so are the errors
(:math:`[\Delta F_\nu] = \mathrm{Jy}`). To perform the line fitting
necessary to determine the spectral index the wavelength at which each
flux density measurement is taken is also needed and must be provided in
microns (:math:`[\lambda] = \mu\mathrm{m}`).

.. autoclass:: cls.classify.star
    :members:
    :private-members: