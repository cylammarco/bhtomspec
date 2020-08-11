.. BHTOM-SPEC documentation master file, created by
   sphinx-quickstart on Thu Jan  9 11:38:33 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

How to Use The BHTOM-SPEC Documentation
=======================================

This provides the spectroscopic data reduction for BHTOM. At the time of writing, it is planned for implementing the longslit spectral reduction of Gemini/GMOS, LT/SPRAT and LCO/Flloyds.


Usage
=====

To perform a spectral extraction from SPRAT data, user only needs to type:

.. code-block:: shell

    python sprat_extraction.py extraction_config.yaml

While for GMOS, the images have to be flattened first, before performing the spectral extraction:

.. code-block:: shell

    python gmos_reduction.py reduction_config.yaml
    python gmos_extraction.py extraction_config.yaml

Some the Configuration Files section for more information.

User Guide
==========

.. toctree::
   :maxdepth: 2
   :caption: Installation guide

   installation

.. toctree::
   :maxdepth: 2
   :caption: Configuration Files

   config/reduction_config
   config/extraction_config

.. toctree::
   :maxdepth: 2
   :caption: Instruments

   examples/ltsprat
   examples/gmosls

Indices and tables
==================

* :ref:`search`
