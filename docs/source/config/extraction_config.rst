.. _extraction_config:

Extraction Configuration
========================

This applies to both GMOS and SPRAT.

Extraction configuration (YAML)
------------------------------------
Once the flattening is performed, the longslit spectral extraction is like any extractions from a generic spectrograph. The format of the parameter files for GMOS and SPRAT are the same. Most of the parameters can be kept constant, respective to the instrument. The file paths, which are the first few lines of the file, have to be updated for each extraction. The output data can be exported as FITS/CSV, static png/jpg/pdf image, an interactive plotly iframe, or just the JSON file for the plotly graph. The following shows the example YAML of the GMOS extraction.

.. code-block:: yaml
   :linenos:

   folder_path: "."
   output_folder_path: "20181115_gaia18cnz/output"

   # Input
   science_file:
     # - they have to be in a list
     - "20181115_gaia18cnz/output/N20181115S0097_flattened.fits"
     - "20181115_gaia18cnz/output/N20181115S0098_flattened.fits"
     - "20181115_gaia18cnz/output/N20181115S0101_flattened.fits"
     - "20181115_gaia18cnz/output/N20181115S0102_flattened.fits"
   science_arc_file: "20181115_gaia18cnz/output/N20181115S0215_flattened.fits"

   standard_file: "20180811_g191b2b/output/N20180811S0148_flattened.fits"
   standard_arc_file: "20180811_g191b2b/output/N20180813S0115_flattened.fits"

   # If sensitivity curve is provided, standard star extraction will be skipped
   sensitivity_file:  "20180811_g191b2b/output/sensitivity_itp.npy"

   # Output - they have to be lists
   science_output_file:
     # - they have to be in a list
     - "20181115_gaia18cnz/output/N20181115S0097"
     - "20181115_gaia18cnz/output/N20181115S0098"
     - "20181115_gaia18cnz/output/N20181115S0101"
     - "20181115_gaia18cnz/output/N20181115S0102"
   standard_output_file: "20180811_g191b2b/output/N20180811S0148"

   science_output_type:
     - "fits"
     - "csv"
     - "png"
     - "iframe"
   standard_output_type:
     - "fits"
     - "csv"
     - "png"
     - "iframe"

   # Output
   savefits: True
   fits_output: 'flux_resampled+wavecal+flux+adu+adu_resampled'
   fits_filename: 'reduced'
   fits_stype: 'science'
   fits_individual: False
   fits_overwrite: False
   savecsv: True
   csv_output: 'flux_resampled+wavecal+flux+adu+adu_resampled'
   csv_filename: 'reduced'
   csv_stype: 'science'
   csv_individual: False
   csv_overwrite: False

   # Only needed if images/iframe are being exported
   width: 1920
   height: 1080

It is then followed by all the extraction parameters. At the time of writing, the readnoise, gain, seeing and exptime are not read automatically from the header. In order to have a correct absolute flux calibration, the exposure time has to be filled in manually. The gain is constant for the instrument over time, so it does not need adusting unless the observatory has changed the instrument setting. While the seeing varies with time, it does not affect the extraction quality, it can slightly slow down the speed in aperture tracing, which is hardly a bottleneck process of the data extraction to begin with, so it can be left as it is. The same applies to the paramteers for the standard.

.. code-block:: yaml
   :linenos:

   # Science traget

   # TwoDSpec parameters
   science_saxis: 1
   science_spatial_mask: ''
   science_spec_mask: ''
   science_flip: True
   science_cosmicray: True
   science_cosmicray_sigma: 5.
   science_readnoise: 0.
   science_gain: 1.64
   science_seeing: 1.2
   science_exptime: 300
   science_silence: False

Below is the rest of the configuration file.

.. code-block:: yaml
   :linenos:

   # Aperture Tracing
   science_aptrace_nspec: 1
   science_aptrace_nwindow: 25
   science_aptrace_spec_sep: 5
   science_aptrace_resample_factor: 10
   science_aptrace_rescale: False
   science_aptrace_scaling_min: 0.995
   science_aptrace_scaling_max: 1.005
   science_aptrace_scaling_step: 0.001
   science_aptrace_percentile: 5
   science_aptrace_tol: 3
   science_aptrace_polydeg: 3
   science_aptrace_ap_faint: 10
   science_aptrace_display: False
   science_aptrace_renderer: "default"
   science_aptrace_jsonstring: False
   science_aptrace_iframe: False
   science_aptrace_open_iframe: False
   # Aperture Extraction
   science_apextract_apwidth: 10
   science_apextract_skysep: 5
   science_apextract_skywidth: 5
   science_apextract_skydeg: 1
   science_apextract_optimal: True
   science_apextract_display: False
   science_apextract_renderer: "default"
   science_apextract_jsonstring: False
   science_apextract_iframe: False
   science_apextract_open_iframe: False
   # WavelengthCalibration
   science_silence: False
   science_pixel_list: ~
   # Provide polyfit coefficients
   science_polyfit: False
   science_polyfit_coeff: [0., 0., 0., 0., 0.]
   science_polyfit_type: ['poly']
   # Extract Arc spectrum
   science_arc_spec_display: False
   science_arc_spec_jsonstring: False
   science_arc_spec_renderer: 'default'
   science_arc_spec_iframe: False
   science_arc_spec_open_iframe: False
   # Find arc lines
   science_findarc_background: 1000.
   science_findarc_percentile: ~
   science_findarc_prominence: 0.
   science_findarc_distance: 5.
   science_findarc_refine: True
   science_findarc_refine_window_width: 5
   science_findarc_display: False
   science_findarc_jsonstring: False
   science_findarc_renderer: 'default'
   science_findarc_iframe: False
   science_findarc_open_iframe: False
   # Calibrator parameters
   science_wavecal_min_wavelength: 5000
   science_wavecal_max_wavelength: 9500
   science_wavecal_plotting_library: 'plotly'
   science_wavecal_log_level: 'info'
   # Calibrator Fit constraints
   science_constraints_num_slopes: 10000
   science_constraints_range_tolerance: 500
   science_constraints_fit_tolerance: 10.
   science_constraints_polydeg: 4
   science_constraints_candidate_thresh: 15.
   science_constraints_linearity_thresh: 1.5
   science_constraints_ransac_thresh: 3
   science_constraints_num_candidates: 25
   science_constraints_xbins: 200
   science_constraints_ybins: 200
   science_constraints_brute_force: False
   science_constraints_polyfit_type: 'poly'
   science_constraints_spec_id: ~
   # Atlas
   science_atlas_elements: ['Cu', 'Ar']
   science_atlas_min_atlas_wavelength: 0
   science_atlas_max_atlas_wavelength: 15000
   science_atlas_min_intensity: 0
   science_atlas_min_distance: 0
   science_atlas_vacuum: False
   science_atlas_pressure: ~
   science_atlas_temperature: ~
   science_atlas_relative_humidity: ~
   science_atlas_constrain_poly: False
   science_atlas_spec_id: ~
   # User-supplied Atlas, the following OVERRIDES the Atlas config set above
   science_atlas_user_supplied: True
   science_atlas_user_wavelengths: [
                   4703.632, 4728.19041, 4766.19677, 4807.36348, 4849.16386, 4881.22627, 4890.40721, 4906.12088, 4934.58593, 4966.46490,
                   5018.56194, 5063.44827, 5163.723, 5189.191, 5497.401,
                   5560.246, 5608.290, 5913.723,
                   6754.698, 6873.185, 6967.352,
                   7032.190, 7069.167, 7149.012, 7274.940, 7386.014,
                   7505.935, 7516.721, 7637.208, 7725.887, 7893.246, 7950.362,
                   8105.921, 8117.542, 8266.794, 8410.521, 8426.963,
                   8523.783, 8670.325,
                   9125.471, 9197.161, 9227.03, 9356.787,
                   9660.435, 9787.186
               ]
   science_atlas_user_elements: ['CuAr']
   science_atlas_user_vacuum: True
   science_atlas_user_pressure: ~
   science_atlas_user_temperature: ~
   science_atlas_user_relative_humidity: ~
   science_atlas_user_constrain_poly: False
   science_atlas_user_spec_id: ~

   # Fit
   science_fit_sample_size: 5
   science_fit_top_n: 10
   science_fit_max_tries: 20000
   science_fit_progress: True
   science_fit_coeff: ~
   science_fit_linear: True
   science_fit_weighted: True
   science_fit_filter_close: False
   science_fit_display: False
   science_fit_savefig: False
   science_fit_filename: ~
   science_fit_spec_id: ~
   # Refine Fit 1st pass
   science_refinefit1: False
   science_refinefit1_polyfit_coeff: ~
   science_refinefit1_n_delta: 2
   science_refinefit1_refine: True
   science_refinefit1_tolerance: 10.
   science_refinefit1_method: 'Nelder-Mead'
   science_refinefit1_convergence: 1.0e-6
   science_refinefit1_robust_refit: True
   science_refinefit1_polydeg: ~
   science_refinefit1_display: False
   science_refinefit1_savefig: False
   science_refinefit1_filename: ~
   science_refinefit1_spec_id: ~
   # Refine Fit 2nd pass
   science_refinefit2: False
   science_refinefit2_polyfit_coeff: ~
   science_refinefit2_n_delta: ~
   science_refinefit2_refine: True
   science_refinefit2_tolerance: 5.
   science_refinefit2_method: 'Nelder-Mead'
   science_refinefit2_convergence: 1.0e-6
   science_refinefit2_robust_refit: True
   science_refinefit2_polydeg: ~
   science_refinefit2_display: False
   science_refinefit2_savefig: False
   science_refinefit2_filename: ~
   science_refinefit2_spec_id: ~

   # Standard traget

   # TwoDSpec parameters
   standard_saxis: 1
   standard_spatial_mask: ''
   standard_spec_mask: ''
   standard_flip: True
   standard_cosmicray: True
   standard_cosmicray_sigma: 5.
   standard_readnoise: ~
   standard_gain: ~
   standard_seeing: ~
   standard_exptime: ~
   standard_silence: False
   # TwoDSpec header keywords (NOT used if the values are provided above)
   standard_readnoise_keyword: ~
   standard_gain_keyword: ~
   standard_seeing_keyword: ~
   standard_exptime_keyword: ~
   # Aperture Tracing
   standard_aptrace_nspec: 1
   standard_aptrace_nwindow: 25
   standard_aptrace_spec_sep: 5
   standard_aptrace_resample_factor: 10
   standard_aptrace_rescale: False
   standard_aptrace_scaling_min: 0.995
   standard_aptrace_scaling_max: 1.005
   standard_aptrace_scaling_step: 0.001
   standard_aptrace_percentile: 5
   standard_aptrace_tol: 3
   standard_aptrace_polydeg: 3
   standard_aptrace_ap_faint: 10
   standard_aptrace_display: False
   standard_aptrace_renderer: "default"
   standard_aptrace_jsonstring: False
   standard_aptrace_iframe: False
   standard_aptrace_open_iframe: False
   # Aperture Extraction
   standard_apextract_apwidth: 15
   standard_apextract_skysep: 5
   standard_apextract_skywidth: 5
   standard_apextract_skydeg: 1
   standard_apextract_optimal: True
   standard_apextract_display: False
   standard_apextract_renderer: "default"
   standard_apextract_jsonstring: False
   standard_apextract_iframe: False
   standard_apextract_open_iframe: False
   # WavelengthCalibration
   standard_silence: False
   standard_pixel_list: ~
   # Provide polyfit coefficients
   standard_polyfit: False
   standard_polyfit_coeff: [0., 0., 0., 0., 0.]
   standard_polyfit_type: ['poly']
   # Extract Arc spectrum
   standard_arc_spec_display: False
   standard_arc_spec_jsonstring: False
   standard_arc_spec_renderer: 'default'
   standard_arc_spec_iframe: False
   standard_arc_spec_open_iframe: False
   # Find arc lines
   standard_findarc_background: 1000.
   standard_findarc_percentile: ~
   standard_findarc_prominence: 0.
   standard_findarc_distance: 5.
   standard_findarc_refine: True
   standard_findarc_refine_window_width: 5
   standard_findarc_display: False
   standard_findarc_jsonstring: False
   standard_findarc_renderer: 'default'
   standard_findarc_iframe: False
   standard_findarc_open_iframe: False
   # Calibrator parameters
   standard_wavecal_min_wavelength: 5000
   standard_wavecal_max_wavelength: 9500
   standard_wavecal_plotting_library: 'plotly'
   standard_wavecal_log_level: 'info'
   # Calibrator Fit constraints
   standard_constraints_num_slopes: 10000
   standard_constraints_range_tolerance: 500
   standard_constraints_fit_tolerance: 10.
   standard_constraints_polydeg: 4
   standard_constraints_candidate_thresh: 15.
   standard_constraints_linearity_thresh: 1.5
   standard_constraints_ransac_thresh: 3
   standard_constraints_num_candidates: 25
   standard_constraints_xbins: 200
   standard_constraints_ybins: 200
   standard_constraints_brute_force: False
   standard_constraints_polyfit_type: 'poly'
   standard_constraints_spec_id: ~
   # Atlas
   standard_atlas_elements: ['Cu', 'Ar']
   standard_atlas_min_atlas_wavelength: 0
   standard_atlas_max_atlas_wavelength: 15000
   standard_atlas_min_intensity: 0
   standard_atlas_min_distance: 0
   standard_atlas_vacuum: False
   standard_atlas_pressure: ~
   standard_atlas_temperature: ~
   standard_atlas_relative_humidity: ~
   standard_atlas_constrain_poly: False
   standard_atlas_spec_id: ~
   # User-supplied Atlas, the following OVERRIDES the Atlas config set above
   standard_atlas_user_supplied: True
   standard_atlas_user_wavelengths: [
                   4703.632, 4728.19041, 4766.19677, 4807.36348, 4849.16386, 4881.22627, 4890.40721, 4906.12088, 4934.58593, 4966.46490,
                   5018.56194, 5063.44827, 5163.723, 5189.191, 5497.401,
                   5560.246, 5608.290, 5913.723,
                   6754.698, 6873.185, 6967.352,
                   7032.190, 7069.167, 7149.012, 7274.940, 7386.014,
                   7505.935, 7516.721, 7637.208, 7725.887, 7893.246, 7950.362,
                   8105.921, 8117.542, 8266.794, 8410.521, 8426.963,
                   8523.783, 8670.325,
                   9125.471, 9197.161, 9227.03, 9356.787,
                   9660.435, 9787.186
               ]
   standard_atlas_user_elements: ['CuAr']
   standard_atlas_user_vacuum: True
   standard_atlas_user_pressure: ~
   standard_atlas_user_temperature: ~
   standard_atlas_user_relative_humidity: ~
   standard_atlas_user_constrain_poly: False
   standard_atlas_user_spec_id: ~
   # Fit
   standard_fit_sample_size: 5
   standard_fit_top_n: 20
   standard_fit_max_tries: 20000
   standard_fit_progress: True
   standard_fit_coeff: ~
   standard_fit_linear: True
   standard_fit_weighted: True
   standard_fit_filter_close: False
   standard_fit_display: False
   standard_fit_savefig: False
   standard_fit_filename: ~
   standard_fit_spec_id: ~
   # Refine Fit 1st pass
   standard_refinefit1: True
   standard_refinefit1_polyfit_coeff: ~
   standard_refinefit1_n_delta: 2
   standard_refinefit1_refine: True
   standard_refinefit1_tolerance: 10.
   standard_refinefit1_method: 'Nelder-Mead'
   standard_refinefit1_convergence: 1.0e-6
   standard_refinefit1_robust_refit: True
   standard_refinefit1_polydeg: ~
   standard_refinefit1_display: False
   standard_refinefit1_savefig: False
   standard_refinefit1_filename: ~
   standard_refinefit1_spec_id: ~
   # Refine Fit 2nd pass
   standard_refinefit2: True
   standard_refinefit2_polyfit_coeff: ~
   standard_refinefit2_n_delta: ~
   standard_refinefit2_refine: True
   standard_refinefit2_tolerance: 5.
   standard_refinefit2_method: 'Nelder-Mead'
   standard_refinefit2_convergence: 1.0e-6
   standard_refinefit2_robust_refit: True
   standard_refinefit2_polydeg: ~
   standard_refinefit2_display: False
   standard_refinefit2_savefig: False
   standard_refinefit2_filename: ~
   standard_refinefit2_spec_id: ~

   # Choose Flux Calibration
   fluxcal_target: 'LTT7987'
   fluxcal_library: 'esoxshooter'
   fluxcal_ftype: 'flux'
   fluxcal_cutoff: 0.4
   fluxcal_display: True
   fluxcal_renderer: 'default'
   fluxcal_jsonstring: False
   fluxcal_iframe: True
   fluxcal_open_iframe: True
   # Compute sensitivity curve
   sensecurve_kind: 3
   sensecurve_smooth: False
   sensecurve_slength: 5
   sensecurve_sorder: 3
   sensecurve_mask_range:
     - [6860, 6960]
     - [7150, 7410]
     - [7580, 7720]
   sensecurve_display: False
   sensecurve_renderer: 'default'
   sensecurve_jsonstring: False
   sensecurve_iframe: False
   sensecurve_open_iframe: False
