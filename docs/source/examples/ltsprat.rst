.. _ltsprat:

********
LT SPRAT
********

The Liverpool Telescope distributes the `SPRAT Level 2 (L2) data <http://telescope.livjm.ac.uk/TelInst/Inst/SPRAT/#filenames>`_ in multi-extension FITS files containing a primary image array for the Level 1 (L1) reduced CCD frame and up to five FITS extensions of derived products. The executive summary of the data product is listed in the following table, full description can be referred to the same link above.

.. list-table:: Multi-extention FITS data format
    :widths: 10 12 60
    :header-rows: 1
    :stub-columns: 1

    * - Index
      - EXTNAME
      - Description

    * - 0
      - L1_IMAGE
      - Fieldflattened long-slit image.

    * - 1
      - LSS_NONSS
      - Wavelength calibrated and trimmed long-slit image.

    * - 2
      - SPEC_NONSS
      - Tophat extracted 1D spectrum.

    * - 3
      - SPEC_SS
      - Tophat extracted and sky subtracted 1D spectrum. Always exists if SPEC_NONSS exists, but it may not be picking the right regions for sky subtraction (e.g. nebulosity from the resolved host galaxy).

    * - 4
      - NORMFLUX
      - Detector response corrected SPEC_SS spectrum. Flux is normalised to unity over the region 5000-6000 A. Always exists if SPEC_SS exists.

    * - 5
      - FLUX
      - Absolute flux calibrated spectrum in unit of erg / s / cm^2 / A.

The provided reduced data only guarantee to contain index 0 and 1. The native pipeline may fail to produce SPEC_NONSS in crowded field; the flux may not be calibrated if the photometry from the acquisition image is not reliable (`see more here <https://github.com/LivTel/sprat_l2_pipeline>`_).

With only the Science FITS file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the data comes with 6 HDUs, the spectrum can be re-extracted directly from the LSS_NONSS image, and then re-appply the absolute flux calibration by extracting the sensitivity function by dividing [SPEC_SS] by [FLUX]. The sensitivity function the native pipeline applied is (ADU / s) /(erg / s / cm^2 / A). This can improve the signal-to-noise ratio (SNR) by up to ~30% coming purely from the difference between optimal extraction and tophat etraction.

If the data comes with 5 HDUs, the same procedure applies except the absolute calibration will be done by using a pre-saved sensitivity function adjusted to the unit of (erg / s / cm^2 / A) / s. This will only give the estimate of the absolute flux. The spectral shape will be the same as the NORMFLUX data from their native pipeline.

If the data comes with less than 5 HDUs, the sensitivity function will be taken from the pre-saved data, which can be out of date but should give a good general shape of the spectrum.

With only the Science and Standard FITS file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is not common for the LT programmes to include standard frames. However, if one is provided, it can be used to compute a more representatibe response function. The procedures are the same as above in terms of spectral extraction. The sensitivity will be computed by providing the name of the standard star in the configuration file.

With arc files
^^^^^^^^^^^^^^
The arc file is taken every time after an observing group/block. It is extremely rare that the automated LT pipeline fail to the perform wavelength calibration. Using a different way to compute the wavelength solution with more lines with RASCAL, a recalibration may give a marginal improvment in the calibration (less than the resolution in any case). We recommend providing an arc file if you wish to use the native pixel scale such that the resampling done from L1_IMAGE to LSS_NONSS will not lead to any loss of information, regardless of how small they can be.

Data Processing
^^^^^^^^^^^^^^^
The above processes will be assignd as one of the following five scenarios when handled by this pipeline:

1. science frame only (Science L1 Image)

   i. Improving the S/N ratio by using optimal extraction (Science)

2. science + science arc frames only (Science LSS_NONSS)

   i. Improving the S/N ratio by using optimal extraction (Science)
   ii. Recomputing the wavelength solution (Science)

3. science + standard frames only (Science L1 Image + Standard L1 Image)

   i. Improving the S/N ratio by using optimal extraction (Science + Standard)
   ii. Recompute the sensitivity response function

4. scince + science arc + standard frames (Science L1 Image + Standard LSS_NONSS)

   i. Improving the S/N ratio by using optimal extraction (Science + Standard)
   ii. Recomputing the wavelength solution (Science)
   iii. Recompute the sensitivity response function

5. scince + science arc + standard + standard arc frames (Science L1 Image + Standard L1 Image)

   i. Improving the S/N ratio by using optimal extraction (Science + Standard)
   ii. Recomputing the wavelength solution (Science + Standard)
   iii. Recompute the sensitivity response function

* science + standard + standard arc will be treated as (3)
* science + standard arc will be treated as (1)
* If sensitivity file is provided, they will be used and the reprocessing will ignore the standard and standard arc files, i.e. (3) will reduced to (1), (4) & (5) will reduced to (2).
