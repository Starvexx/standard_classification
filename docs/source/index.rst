.. YSO standard classification documentation master file, created by
   sphinx-quickstart on Wed Jan 10 12:48:07 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to YSO standard classification's documentation!
=======================================================

This Python package is used to classify young stellar objects using the
standard classification scheme originally devised by Lada et al. (1984).

It determines the infrared spectral index :math:`\alpha_\mathrm{IR}` in the range
from :math:`2\,\mu\mathrm{m}` to :math:`20\,\mu\mathrm{m}`, estimating the slope
of the spectral energy distribution. Based on the measured slope, the
YSOs are divided into 5 distinct classes (0/I, flat spectrum, II,
III with debris disk, III no disk/MS star).


.. note::
   This project is under active development.

Contents
--------

.. toctree::
   api/star