import os
import sys
import warnings

from aspired import image_reduction
from aspired import spectral_reduction
from astropy.io import fits
import numpy as np
import yaml

from .gmos_fieldflattening import create_pixel_array

base_path = os.path.abspath(os.path.dirname(__file__))

# Get config file
try:
    params_path = os.path.join(base_path, sys.argv[1])
except:
    params_path = os.path.join(base_path, 'extraction_config.yaml')

if not os.path.isabs(params_path):
    params_path = os.path.abspath(params_path)

with open(params_path, 'r') as stream:
    params = yaml.safe_load(stream)
    print('Reading parameters from ' + params_path)

# Get the working directory and output directory
params_folder_path = os.path.dirname(params_path)
folder_path = os.path.join(params_folder_path,
                           os.path.dirname(params['folder_path']))
output_folder_path = os.path.join(
    params_folder_path, os.path.dirname(params['output_folder_path']))

if not os.path.isabs(folder_path):
    folder_path = os.path.abspath(folder_path)

if output_folder_path is None:
    output_folder_path = folder_path

if not os.path.isabs(output_folder_path):
    output_folder_path = os.path.join(folder_path, output_folder_path)

if not os.path.exists(output_folder_path):
    os.mkdir(output_folder_path)

if params['science_file'] is None:
    ValueError('Science frame is not provided. Reduction terminates.')


# Get the filepaths
science_filepath = [
    os.path.join(folder_path, sci_path) for sci_path in params['science_file']
]

if params['science_arc_file'] is not None:
    science_arc_filepath = os.path.join(folder_path,
                                        params['science_arc_file'])
else:
    science_arc_filepath = None


if params['standard_file'] is not None:
    standard_filepath = os.path.join(folder_path, params['standard_file'])
else:
    standard_filepath = None

if params['standard_arc_file'] is not None:
    standard_arc_filepath = os.path.join(folder_path,
                                         params['standard_arc_file'])
else:
    standard_arc_filepath = None

if params['sensitivity_file'] is not None:
    sensitivity_filepath = os.path.join(folder_path,
                                        params['sensitivity_file'])
else:
    sensitivity_filepath = None

# Get the output filepaths and types
science_output_filepath = [
    os.path.join(output_folder_path, sci_path)
    for sci_path in params['science_output_file']
]
science_output_type = params['science_output_type']

try:
    standard_output_filepath = os.path.join(output_folder_path,
                                            params['standard_output_file'])
    standard_output_type = params['standard_output_type']
except:
    pass

'''
The process will be assignd as one of the four scenarios, there must be at
least one science frame.

(1) science frame only
    i. Count as a function of pixel in the dispersion direction

(2) science + science arc frames only
    i. Counts as a function of wavelength

(3a) scince + science arc + standard frames
    i. Flux as a function of wavelength

(3b) scince + science arc + standard + standard arc frames
    i. Flux as a function of wavelength (more accurate than case 3)

(*) science + standard + standard arc will be treated as (3)
(**) science + standard arc will be treated as (2)
(***) science + standard frames only will be treated as (1)

If sensitivity file is provided, they will be used and the reprocessing will
ignore the standard and standard arc files, i.e. the processing of (3) & (4)
will be reduced to (2).

To summarise in a table:

+------+---------+-------------+----------+--------------+
| Case | Science | Science Arc | Standard | Standard Arc |
+------+---------+-------------+----------+--------------+
|  1.  |    v    |      x      |     x    |      x       |
|  2.  |    v    |      v      |     x    |      x       |
| (3.) |    v    |      v      |     v    |      x       |
|  3.  |    v    |      v      |     v    |      v       |
+------+---------+-------------+----------+--------------+
| (1.) |    v    |      x      |     v    |      x       |
| (2.) |    v    |      x      |     x    |      v       |
| (3.) |    v    |      x      |     v    |      v       |
+------+---------+-------------+----------+--------------+

'''

# with science arc
if (science_arc_filepath is not None):
    # with standard
    if (standard_filepath is not None):
        # with standard arc
        if (standard_arc_filepath is not None):
            case = 3
        # without standard arc
        else:
            standard_arc_filepath = science_arc_filepath
            case = 3
    # without standard
    else:
        case = 2
# without science arc
else:
    # with standard
    if (standard_filepath is not None):
        # with standard arc
        if (standard_arc_filepath is not None):
            case = 3
            science_arc_filepath = standard_arc_filepath
        # without standard arc
        else:
            science_filepath.append(standard_filepath)
            science_output_filepath.append(standard_output_filepath)
            science_output_type.append(standard_output_type)
            case = 1
    # without standard
    else:
        # with standard arc
        if (standard_arc_filepath is not None):
            science_arc_filepath = standard_arc_filepath
            case = 2
        # without standard arc
        else:
            case = 1

if "png" in science_output_type:
    save_png = True
else:
    save_png = False

if "jpg" in science_output_type:
    save_jpg = True
else:
    save_jpg = False

if "svg" in science_output_type:
    save_svg = True
else:
    save_svg = False

if "pdf" in science_output_type:
    save_pdf = True
else:
    save_pdf = False

if "json" in science_output_type:
    save_json = True
else:
    save_json = False

if "iframe" in science_output_type:
    save_iframe = True
else:
    save_iframe = False

if case >= 3:
    if "png" in standard_output_type:
        save_png_standard = True
    else:
        save_png_standard = False

    if "jpg" in standard_output_type:
        save_jpg_standard = True
    else:
        save_jpg_standard = False

    if "svg" in standard_output_type:
        save_svg_standard = True
    else:
        save_svg_standard = False

    if "pdf" in standard_output_type:
        save_pdf_standard = True
    else:
        save_pdf_standard = False

    if "json" in standard_output_type:
        save_json_standard = True
    else:
        save_json_standard = False

    if "iframe" in standard_output_type:
        save_iframe_standard = True
    else:
        save_iframe_standard = False

# Twodspec
science_saxis = params['science_saxis']
science_spatial_mask = params['science_spatial_mask']
science_spec_mask = params['science_spec_mask']
science_flip = params['science_flip']
science_cosmicray = params['science_cosmicray']
science_cosmicray_sigma = params['science_cosmicray_sigma']
science_readnoise = params['science_readnoise']
science_gain = params['science_gain']
science_seeing = params['science_seeing']
science_exptime = params['science_exptime']
science_silence = params['science_silence']

# Get the keywords (NOT used if the values are provided above)
science_readnoise_keyword = params['science_readnoise_keyword']
science_gain_keyword = params['science_gain_keyword']
science_seeing_keyword = params['science_seeing_keyword']
science_exptime_keyword = params['science_exptime_keyword']

# Get the aperture tracing parameters
science_aptrace_nspec = params['science_aptrace_nspec']
science_aptrace_nwindow = params['science_aptrace_nwindow']
science_aptrace_spec_sep = params['science_aptrace_spec_sep']
science_aptrace_resample_factor = params['science_aptrace_resample_factor']
science_aptrace_rescale = params['science_aptrace_rescale']
science_aptrace_scaling_min = params['science_aptrace_scaling_min']
science_aptrace_scaling_max = params['science_aptrace_scaling_max']
science_aptrace_scaling_step = params['science_aptrace_scaling_step']
science_aptrace_percentile = params['science_aptrace_percentile']
science_aptrace_tol = params['science_aptrace_tol']
science_aptrace_polydeg = params['science_aptrace_polydeg']
science_aptrace_ap_faint = params['science_aptrace_ap_faint']
science_aptrace_display = params['science_aptrace_display']
science_aptrace_renderer = params['science_aptrace_renderer']
science_aptrace_jsonstring = params['science_aptrace_jsonstring']
science_aptrace_iframe = params['science_aptrace_iframe']
science_aptrace_open_iframe = params['science_aptrace_open_iframe']

# Get the aperture extract parameters
science_apextract_apwidth = params['science_apextract_apwidth']
science_apextract_skysep = params['science_apextract_skysep']
science_apextract_skywidth = params['science_apextract_skywidth']
science_apextract_skydeg = params['science_apextract_skydeg']
science_apextract_optimal = params['science_apextract_optimal']
science_apextract_display = params['science_apextract_display']
science_apextract_renderer = params['science_apextract_renderer']
science_apextract_jsonstring = params['science_apextract_jsonstring']
science_apextract_iframe = params['science_apextract_iframe']
science_apextract_open_iframe = params['science_apextract_open_iframe']

# WavelengthCalibration
science_silence = params['science_silence']
science_pixel_list = params['science_pixel_list']

# Extract Arc spectrum
science_arc_spec_display = params['science_arc_spec_display']
science_arc_spec_jsonstring = params['science_arc_spec_jsonstring']
science_arc_spec_renderer = params['science_arc_spec_renderer']
science_arc_spec_iframe = params['science_arc_spec_iframe']
science_arc_spec_open_iframe = params['science_arc_spec_open_iframe']

# Find arc lines
science_findarc_background = params['science_findarc_background']
science_findarc_percentile = params['science_findarc_percentile']
science_findarc_prominence = params['science_findarc_prominence']
science_findarc_distance = params['science_findarc_distance']
science_findarc_refine = params['science_findarc_refine']
science_findarc_refine_window_width = params[
    'science_findarc_refine_window_width']
science_findarc_display = params['science_findarc_display']
science_findarc_jsonstring = params['science_findarc_jsonstring']
science_findarc_renderer = params['science_findarc_renderer']
science_findarc_iframe = params['science_findarc_iframe']
science_findarc_open_iframe = params['science_findarc_open_iframe']

# Calibrator parameters
science_wavecal_min_wavelength = params['science_wavecal_min_wavelength']
science_wavecal_max_wavelength = params['science_wavecal_max_wavelength']
science_wavecal_plotting_library = params['science_wavecal_plotting_library']
science_wavecal_log_level = params['science_wavecal_log_level']

# Provide polyfit coefficients
science_wavecal_polyfit = params['science_polyfit']
science_wavecal_polyfit_coeff = params['science_polyfit_coeff']
science_wavecal_polyfit_type = params['science_polyfit_type']

# Calibrator Fit constraints
science_constraints_num_slopes = params['science_constraints_num_slopes']
science_constraints_range_tolerance = params[
    'science_constraints_range_tolerance']
science_constraints_fit_tolerance = params['science_constraints_fit_tolerance']
science_constraints_polydeg = params['science_constraints_polydeg']
science_constraints_candidate_thresh = params[
    'science_constraints_candidate_thresh']
science_constraints_linearity_thresh = params[
    'science_constraints_linearity_thresh']
science_constraints_ransac_thresh = params['science_constraints_ransac_thresh']
science_constraints_num_candidates = params[
    'science_constraints_num_candidates']
science_constraints_xbins = params['science_constraints_xbins']
science_constraints_ybins = params['science_constraints_ybins']
science_constraints_brute_force = params['science_constraints_brute_force']
science_constraints_polyfit_type = params['science_constraints_polyfit_type']
science_constraints_spec_id = params['science_constraints_spec_id']

# Atlas
science_atlas_elements = params['science_atlas_elements']
science_atlas_min_atlas_wavelength = params[
    'science_atlas_min_atlas_wavelength']
science_atlas_max_atlas_wavelength = params[
    'science_atlas_max_atlas_wavelength']
science_atlas_min_intensity = params['science_atlas_min_intensity']
science_atlas_min_distance = params['science_atlas_min_distance']
science_atlas_vacuum = params['science_atlas_vacuum']
science_atlas_pressure = params['science_atlas_pressure']
science_atlas_temperature = params['science_atlas_temperature']
science_atlas_relative_humidity = params['science_atlas_relative_humidity']
science_atlas_constrain_poly = params['science_atlas_constrain_poly']
science_atlas_spec_id = params['science_atlas_spec_id']

# Atlas - User supplied
science_atlas_user_supplied = params['science_atlas_user_supplied']
science_atlas_user_wavelengths = params['science_atlas_user_wavelengths']
science_atlas_user_elements = params['science_atlas_user_elements']
science_atlas_user_vacuum = params['science_atlas_user_vacuum']
science_atlas_user_pressure = params['science_atlas_user_pressure']
science_atlas_user_temperature = params['science_atlas_user_temperature']
science_atlas_user_relative_humidity = params[
    'science_atlas_user_relative_humidity']
science_atlas_user_constrain_poly = params['science_atlas_user_constrain_poly']
science_atlas_user_spec_id = params['science_atlas_user_spec_id']

# Fit
science_fit_sample_size = params['science_fit_sample_size']
science_fit_top_n = params['science_fit_top_n']
science_fit_max_tries = params['science_fit_max_tries']
science_fit_progress = params['science_fit_progress']
science_fit_coeff = params['science_fit_coeff']
science_fit_linear = params['science_fit_linear']
science_fit_weighted = params['science_fit_weighted']
science_fit_filter_close = params['science_fit_filter_close']
science_fit_display = params['science_fit_display']
science_fit_savefig = params['science_fit_savefig']
science_fit_filename = params['science_fit_filename']
science_fit_spec_id = params['science_fit_spec_id']

# Refine Fit 1st pass
science_refinefit1 = params['science_refinefit1']
science_refinefit1_polyfit_coeff = params['science_refinefit1_polyfit_coeff']
science_refinefit1_n_delta = params['science_refinefit1_n_delta']
science_refinefit1_refine = params['science_refinefit1_refine']
science_refinefit1_tolerance = params['science_refinefit1_tolerance']
science_refinefit1_method = params['science_refinefit1_method']
science_refinefit1_convergence = params['science_refinefit1_convergence']
science_refinefit1_robust_refit = params['science_refinefit1_robust_refit']
science_refinefit1_polydeg = params['science_refinefit1_polydeg']
science_refinefit1_display = params['science_refinefit1_display']
science_refinefit1_savefig = params['science_refinefit1_savefig']
science_refinefit1_filename = params['science_refinefit1_filename']
science_refinefit1_spec_id = params['science_refinefit1_spec_id']

# Refine Fit 2nd pass
science_refinefit2 = params['science_refinefit2']
science_refinefit2_polyfit_coeff = params['science_refinefit2_polyfit_coeff']
science_refinefit2_n_delta = params['science_refinefit2_n_delta']
science_refinefit2_refine = params['science_refinefit2_refine']
science_refinefit2_tolerance = params['science_refinefit2_tolerance']
science_refinefit2_method = params['science_refinefit2_method']
science_refinefit2_convergence = params['science_refinefit2_convergence']
science_refinefit2_robust_refit = params['science_refinefit2_robust_refit']
science_refinefit2_polydeg = params['science_refinefit2_polydeg']
science_refinefit2_display = params['science_refinefit2_display']
science_refinefit2_savefig = params['science_refinefit2_savefig']
science_refinefit2_filename = params['science_refinefit2_filename']
science_refinefit2_spec_id = params['science_refinefit2_spec_id']

# Standard

# Twodspec
standard_saxis = params['standard_saxis']
standard_spatial_mask = params['standard_spatial_mask']
standard_spec_mask = params['standard_spec_mask']
standard_flip = params['standard_flip']
standard_cosmicray = params['standard_cosmicray']
standard_cosmicray_sigma = params['standard_cosmicray_sigma']
standard_readnoise = params['standard_readnoise']
standard_gain = params['standard_gain']
standard_seeing = params['standard_seeing']
standard_exptime = params['standard_exptime']
standard_silence = params['standard_silence']

# Get the keywords (NOT used if the values are provided above)
standard_readnoise_keyword = params['standard_readnoise_keyword']
standard_gain_keyword = params['standard_gain_keyword']
standard_seeing_keyword = params['standard_seeing_keyword']
standard_exptime_keyword = params['standard_exptime_keyword']

# Get the aperture tracing parameters
standard_aptrace_nspec = params['standard_aptrace_nspec']
standard_aptrace_nwindow = params['standard_aptrace_nwindow']
standard_aptrace_spec_sep = params['standard_aptrace_spec_sep']
standard_aptrace_resample_factor = params['standard_aptrace_resample_factor']
standard_aptrace_rescale = params['standard_aptrace_rescale']
standard_aptrace_scaling_min = params['standard_aptrace_scaling_min']
standard_aptrace_scaling_max = params['standard_aptrace_scaling_max']
standard_aptrace_scaling_step = params['standard_aptrace_scaling_step']
standard_aptrace_percentile = params['standard_aptrace_percentile']
standard_aptrace_tol = params['standard_aptrace_tol']
standard_aptrace_polydeg = params['standard_aptrace_polydeg']
standard_aptrace_ap_faint = params['standard_aptrace_ap_faint']
standard_aptrace_display = params['standard_aptrace_display']
standard_aptrace_renderer = params['standard_aptrace_renderer']
standard_aptrace_jsonstring = params['standard_aptrace_jsonstring']
standard_aptrace_iframe = params['standard_aptrace_iframe']
standard_aptrace_open_iframe = params['standard_aptrace_open_iframe']

# Get the aperture extract parameters
standard_apextract_apwidth = params['standard_apextract_apwidth']
standard_apextract_skysep = params['standard_apextract_skysep']
standard_apextract_skywidth = params['standard_apextract_skywidth']
standard_apextract_skydeg = params['standard_apextract_skydeg']
standard_apextract_optimal = params['standard_apextract_optimal']
standard_apextract_display = params['standard_apextract_display']
standard_apextract_renderer = params['standard_apextract_renderer']
standard_apextract_jsonstring = params['standard_apextract_jsonstring']
standard_apextract_iframe = params['standard_apextract_iframe']
standard_apextract_open_iframe = params['standard_apextract_open_iframe']

# WavelengthCalibration
standard_silence = params['standard_silence']
standard_pixel_list = params['standard_pixel_list']

# Extract Arc spectrum
standard_arc_spec_display = params['standard_arc_spec_display']
standard_arc_spec_jsonstring = params['standard_arc_spec_jsonstring']
standard_arc_spec_renderer = params['standard_arc_spec_renderer']
standard_arc_spec_iframe = params['standard_arc_spec_iframe']
standard_arc_spec_open_iframe = params['standard_arc_spec_open_iframe']

# Find arc lines
standard_findarc_background = params['standard_findarc_background']
standard_findarc_percentile = params['standard_findarc_percentile']
standard_findarc_prominence = params['standard_findarc_prominence']
standard_findarc_distance = params['standard_findarc_distance']
standard_findarc_refine = params['standard_findarc_refine']
standard_findarc_refine_window_width = params[
    'standard_findarc_refine_window_width']
standard_findarc_display = params['standard_findarc_display']
standard_findarc_jsonstring = params['standard_findarc_jsonstring']
standard_findarc_renderer = params['standard_findarc_renderer']
standard_findarc_iframe = params['standard_findarc_iframe']
standard_findarc_open_iframe = params['standard_findarc_open_iframe']

# Calibrator parameters
standard_wavecal_min_wavelength = params['standard_wavecal_min_wavelength']
standard_wavecal_max_wavelength = params['standard_wavecal_max_wavelength']
standard_wavecal_plotting_library = params['standard_wavecal_plotting_library']
standard_wavecal_log_level = params['standard_wavecal_log_level']

# Provide polyfit coefficients
standard_wavecal_polyfit = params['standard_polyfit']
standard_wavecal_polyfit_coeff = params['standard_polyfit_coeff']
standard_wavecal_polyfit_type = params['standard_polyfit_type']

# Calibrator Fit constraints
standard_constraints_num_slopes = params['standard_constraints_num_slopes']
standard_constraints_range_tolerance = params[
    'standard_constraints_range_tolerance']
standard_constraints_fit_tolerance = params[
    'standard_constraints_fit_tolerance']
standard_constraints_polydeg = params['standard_constraints_polydeg']
standard_constraints_candidate_thresh = params[
    'standard_constraints_candidate_thresh']
standard_constraints_linearity_thresh = params[
    'standard_constraints_linearity_thresh']
standard_constraints_ransac_thresh = params[
    'standard_constraints_ransac_thresh']
standard_constraints_num_candidates = params[
    'standard_constraints_num_candidates']
standard_constraints_xbins = params['standard_constraints_xbins']
standard_constraints_ybins = params['standard_constraints_ybins']
standard_constraints_brute_force = params['standard_constraints_brute_force']
standard_constraints_polyfit_type = params['standard_constraints_polyfit_type']
standard_constraints_spec_id = params['standard_constraints_spec_id']

# Atlas
standard_atlas_elements = params['standard_atlas_elements']
standard_atlas_min_atlas_wavelength = params[
    'standard_atlas_min_atlas_wavelength']
standard_atlas_max_atlas_wavelength = params[
    'standard_atlas_max_atlas_wavelength']
standard_atlas_min_intensity = params['standard_atlas_min_intensity']
standard_atlas_min_distance = params['standard_atlas_min_distance']
standard_atlas_vacuum = params['standard_atlas_vacuum']
standard_atlas_pressure = params['standard_atlas_pressure']
standard_atlas_temperature = params['standard_atlas_temperature']
standard_atlas_relative_humidity = params['standard_atlas_relative_humidity']
standard_atlas_constrain_poly = params['standard_atlas_constrain_poly']
standard_atlas_spec_id = params['standard_atlas_spec_id']

# Atlas - User supplied
standard_atlas_user_supplied = params['standard_atlas_user_supplied']
standard_atlas_user_wavelengths = params['standard_atlas_user_wavelengths']
standard_atlas_user_elements = params['standard_atlas_user_elements']
standard_atlas_user_vacuum = params['standard_atlas_user_vacuum']
standard_atlas_user_pressure = params['standard_atlas_user_pressure']
standard_atlas_user_temperature = params['standard_atlas_user_temperature']
standard_atlas_user_relative_humidity = params[
    'standard_atlas_user_relative_humidity']
standard_atlas_user_constrain_poly = params[
    'standard_atlas_user_constrain_poly']
standard_atlas_user_spec_id = params['standard_atlas_user_spec_id']

# Fit
standard_fit_sample_size = params['standard_fit_sample_size']
standard_fit_top_n = params['standard_fit_top_n']
standard_fit_max_tries = params['standard_fit_max_tries']
standard_fit_progress = params['standard_fit_progress']
standard_fit_coeff = params['standard_fit_coeff']
standard_fit_linear = params['standard_fit_linear']
standard_fit_weighted = params['standard_fit_weighted']
standard_fit_filter_close = params['standard_fit_filter_close']
standard_fit_display = params['standard_fit_display']
standard_fit_savefig = params['standard_fit_savefig']
standard_fit_filename = params['standard_fit_filename']
standard_fit_spec_id = params['standard_fit_spec_id']

# Refine Fit 1st pass
standard_refinefit1 = params['standard_refinefit1']
standard_refinefit1_polyfit_coeff = params['standard_refinefit1_polyfit_coeff']
standard_refinefit1_n_delta = params['standard_refinefit1_n_delta']
standard_refinefit1_refine = params['standard_refinefit1_refine']
standard_refinefit1_tolerance = params['standard_refinefit1_tolerance']
standard_refinefit1_method = params['standard_refinefit1_method']
standard_refinefit1_convergence = params['standard_refinefit1_convergence']
standard_refinefit1_robust_refit = params['standard_refinefit1_robust_refit']
standard_refinefit1_polydeg = params['standard_refinefit1_polydeg']
standard_refinefit1_display = params['standard_refinefit1_display']
standard_refinefit1_savefig = params['standard_refinefit1_savefig']
standard_refinefit1_filename = params['standard_refinefit1_filename']
standard_refinefit1_spec_id = params['standard_refinefit1_spec_id']

# Refine Fit 2nd pass
standard_refinefit2 = params['standard_refinefit2']
standard_refinefit2_polyfit_coeff = params['standard_refinefit2_polyfit_coeff']
standard_refinefit2_n_delta = params['standard_refinefit2_n_delta']
standard_refinefit2_refine = params['standard_refinefit2_refine']
standard_refinefit2_tolerance = params['standard_refinefit2_tolerance']
standard_refinefit2_method = params['standard_refinefit2_method']
standard_refinefit2_convergence = params['standard_refinefit2_convergence']
standard_refinefit2_robust_refit = params['standard_refinefit2_robust_refit']
standard_refinefit2_polydeg = params['standard_refinefit2_polydeg']
standard_refinefit2_display = params['standard_refinefit2_display']
standard_refinefit2_savefig = params['standard_refinefit2_savefig']
standard_refinefit2_filename = params['standard_refinefit2_filename']
standard_refinefit2_spec_id = params['standard_refinefit2_spec_id']

# Flux calibration parameters
fluxcal_target = params['fluxcal_target']
fluxcal_library = params['fluxcal_library']
fluxcal_ftype = params['fluxcal_ftype']
fluxcal_cutoff = params['fluxcal_cutoff']
fluxcal_display = params['fluxcal_display']
fluxcal_renderer = params['fluxcal_renderer']
fluxcal_jsonstring = params['fluxcal_jsonstring']
fluxcal_iframe = params['fluxcal_iframe']
fluxcal_open_iframe = params['fluxcal_open_iframe']

# Sensitivity curve paramters
sensecurve_kind = params['sensecurve_kind']
sensecurve_smooth = params['sensecurve_smooth']
sensecurve_slength = params['sensecurve_slength']
sensecurve_sorder = params['sensecurve_sorder']
sensecurve_mask_range = params['sensecurve_mask_range']
sensecurve_mask_fit_order = params['sensecurve_mask_fit_order']
sensecurve_mask_fit_size = params['sensecurve_mask_fit_size']
sensecurve_display = params['sensecurve_display']
sensecurve_renderer = params['sensecurve_renderer']
sensecurve_jsonstring = params['sensecurve_jsonstring']
sensecurve_iframe = params['sensecurve_iframe']
sensecurve_open_iframe = params['sensecurve_open_iframe']

# Output
savefits = params['savefits']
fits_output = params['fits_output']
fits_filename = params['fits_filename']
fits_stype = params['fits_stype']
fits_individual = params['fits_individual']
fits_overwrite = params['fits_overwrite']
savecsv = params['savecsv']
csv_output = params['csv_output']
csv_filename = params['csv_filename']
csv_stype = params['csv_stype']
csv_individual = params['csv_individual']
csv_overwrite = params['csv_overwrite']
width = params['width']
height = params['height']


n_science = len(science_filepath)

science_twodspec = np.array([None] * n_science, dtype='object')
onedspec = np.array([None] * n_science, dtype='object')

for i in range(n_science):

    # Science frames must exist
    science_fits = fits.open(science_filepath[i])[0]
    science_flattened = science_fits.data
    science_header = science_fits.header
    science_northsouth = science_header['NS']
    science_binx = science_header['BINX']
    science_imported = True

    # Initialise the 1D spec
    onedspec[i] = spectral_reduction.OneDSpec()

    # Dealing with spectral extraction and calibration with ASPIRED, only runs
    # if ASPIRED is imported successfully.
    if science_pixel_list is None:
        science_pixel_list = create_pixel_array(science_northsouth,
                                                science_binx)

    if science_atlas_relative_humidity is None:
        science_atlas_relative_humidity = float(science_header['HUMIDITY'])
    if science_atlas_temperature is None:
        science_atlas_temperature = float(
            science_header['TAMBIENT']) + 273.15
    if science_atlas_pressure is None:
        science_atlas_pressure = float(science_header['PRESSUR2'])

    if science_atlas_user_relative_humidity is None:
        science_atlas_user_relative_humidity = float(
            science_header['HUMIDITY'])
    if science_atlas_user_temperature is None:
        science_atlas_user_temperature = float(
            science_header['TAMBIENT']) + 273.15
    if science_atlas_user_pressure is None:
        science_atlas_user_pressure = float(science_header['PRESSUR2'])

    # Make sure standard_noise, gain, seeing, exptime are None for the
    # above keywords to be used.
    science_twodspec[i] = spectral_reduction.TwoDSpec(
        science_flattened.copy(),
        saxis=science_saxis,
        spatial_mask=science_spatial_mask,
        spec_mask=science_spec_mask,
        flip=science_flip,
        cosmicray=science_cosmicray,
        cosmicray_sigma=science_cosmicray_sigma,
        readnoise=science_readnoise,
        gain=science_gain,
        seeing=science_seeing,
        exptime=science_exptime,
        silence=science_silence)

    science_twodspec[i].ap_trace(
        nspec=science_aptrace_nspec,
        nwindow=science_aptrace_nwindow,
        spec_sep=science_aptrace_spec_sep,
        resample_factor=science_aptrace_resample_factor,
        rescale=science_aptrace_rescale,
        scaling_min=science_aptrace_scaling_min,
        scaling_max=science_aptrace_scaling_max,
        scaling_step=science_aptrace_scaling_step,
        percentile=science_aptrace_percentile,
        tol=science_aptrace_tol,
        polydeg=science_aptrace_polydeg,
        ap_faint=science_aptrace_ap_faint,
        display=science_aptrace_display,
        renderer=science_aptrace_renderer,
        return_jsonstring=science_aptrace_jsonstring,
        save_iframe=science_aptrace_iframe,
        open_iframe=science_aptrace_open_iframe)

    science_twodspec[i].ap_extract(
        apwidth=science_apextract_apwidth,
        skysep=science_apextract_skysep,
        skywidth=science_apextract_skywidth,
        skydeg=science_apextract_skydeg,
        optimal=science_apextract_optimal,
        display=science_apextract_display,
        renderer=science_apextract_renderer,
        return_jsonstring=science_apextract_jsonstring,
        save_iframe=science_apextract_iframe,
        open_iframe=science_apextract_open_iframe)
    onedspec[i].from_twodspec(science_twodspec[i],
                             pixel_list=science_pixel_list,
                             stype='science')

    #################################################################

    # # # ##### #   # ##### #     ##### #   # ##### ##### #   #
    # # # #   # #   # #     #     #     ##  # #       #   #   #
    # # # ##### #   # ##### #     ##### # # # # ###   #   #####
    # # # #   #  # #  #     #     #     #  ## #   #   #   #   #
    ##### #   #   #   ##### ##### ##### #   # #####   #   #   #

    ##### ##### #     ##### ####  ####  ##### ##### ##### ##### #   #
    #     #   # #       #   #   # #   # #   #   #     #   #   # ##  #
    #     ##### #       #   ####  ####  #####   #     #   #   # # # #
    #     #   # #       #   #   # #  #  #   #   #     #   #   # #  ##
    ##### #   # ##### ##### ####  #   # #   #   #   ##### ##### #   #

    #################################################################

    if case >= 2:
        science_arc_fits = fits.open(science_arc_filepath)[0]
        science_arc_flattened = science_arc_fits.data
        science_arc_header = science_arc_fits.header
        science_arc_northsouth = science_arc_header['NS']
        science_arc_binx = science_arc_header['BINX']

        # Check arc binning matches other frames
        # GMOS-South didn't take any CuAr R400 in 1x1 binning since 2018
        if science_arc_binx == science_binx:
            onedspec[i].add_arc(np.flip(science_arc_flattened),
                                stype='science')
            onedspec[i].extract_arc_spec(
                display=science_arc_spec_display,
                return_jsonstring=science_arc_spec_jsonstring,
                renderer=science_arc_spec_renderer,
                save_iframe=science_arc_spec_iframe,
                open_iframe=science_arc_spec_open_iframe,
                stype='science')

        else:
            binx_ratio = science_arc_binx / science_binx
            for j in range(science_aptrace_nspec):
                spec = onedspec.science_twodspec[j]

                len_trace = len(spec.trace)
                trace = np.nanmean(spec.trace) / binx_ratio
                trace_width = np.nanmean(spec.trace_sigma) * 3. / binx_ratio

                arc_trace = np.flip(science_arc_flattened)[
                    max(0, int(trace - trace_width -
                               1)):min(int(trace + trace_width), len_trace), :]
                arc_spec = np.nanmedian(arc_trace, axis=0)
                arc_spec_resampled = signal.resample(arc_spec, len_trace)

                onedspec[i].add_arc_spec(spec_id=j,
                                         arc_spec=arc_spec_resampled)

        onedspec[i].find_arc_lines(
            background=science_findarc_background,
            percentile=science_findarc_percentile,
            prominence=science_findarc_prominence,
            distance=science_findarc_distance,
            refine=science_findarc_refine,
            refine_window_width=science_findarc_refine_window_width,
            display=science_findarc_display,
            return_jsonstring=science_findarc_jsonstring,
            renderer=science_findarc_renderer,
            save_iframe=science_findarc_iframe,
            open_iframe=science_findarc_open_iframe,
            stype='science')

        onedspec[i].initialise_calibrator(
            pixel_list=science_pixel_list,
            min_wavelength=science_wavecal_min_wavelength,
            max_wavelength=science_wavecal_max_wavelength,
            plotting_library=science_wavecal_plotting_library,
            log_level=science_wavecal_log_level,
            stype='science')
        if science_wavecal_polyfit:
            onedspec[i].add_polyfit(
                polyfit_coeff=science_wavecal_polyfit_coeff,
                polyfit_type=science_wavecal_polyfit_type,
                stype='science')
        else:
            onedspec[i].set_fit_constraints(
                num_slopes=science_constraints_num_slopes,
                range_tolerance=science_constraints_range_tolerance,
                fit_tolerance=science_constraints_fit_tolerance,
                polydeg=science_constraints_polydeg,
                candidate_thresh=science_constraints_candidate_thresh,
                linearity_thresh=science_constraints_linearity_thresh,
                ransac_thresh=science_constraints_ransac_thresh,
                num_candidates=science_constraints_num_candidates,
                xbins=science_constraints_xbins,
                ybins=science_constraints_ybins,
                brute_force=science_constraints_brute_force,
                polyfit_type=science_constraints_polyfit_type,
                spec_id=science_constraints_spec_id,
                stype='science')
            if science_atlas_user_supplied:
                if len(science_atlas_user_elements) == 1:
                    elements = [science_atlas_user_elements
                                ] * len(science_atlas_user_wavelengths)
                elif len(science_atlas_user_elements) == len(
                        science_atlas_user_wavelengths):
                    elements = science_atlas_user_elements
                else:
                    raise ValueError('The science_atlas_user_elements should '
                                     'have length 1 or the same length as '
                                     'science_atlas_user_wavelengths')
                onedspec[i].load_user_atlas(
                    elements=elements,
                    wavelengths=science_atlas_user_wavelengths,
                    constrain_poly=science_atlas_user_constrain_poly,
                    vacuum=science_atlas_user_vacuum,
                    pressure=science_atlas_user_pressure,
                    temperature=science_atlas_user_temperature,
                    relative_humidity=science_atlas_user_relative_humidity,
                    stype='science')
            else:
                onedspec[i].add_atlas(
                    elements=science_atlas_elements,
                    min_atlas_wavelength=science_atlas_min_atlas_wavelength,
                    max_atlas_wavelength=science_atlas_max_atlas_wavelength,
                    min_intensity=science_atlas_min_intensity,
                    min_distance=science_atlas_min_distance,
                    vacuum=science_atlas_vacuum,
                    pressure=science_atlas_pressure,
                    temperature=science_atlas_temperature,
                    relative_humidity=science_atlas_relative_humidity,
                    constrain_poly=science_atlas_constrain_poly,
                    spec_id=science_atlas_spec_id,
                    stype='science')

            # wavelength solution for target
            onedspec[i].fit(sample_size=science_fit_sample_size,
                            top_n=science_fit_top_n,
                            max_tries=science_fit_max_tries,
                            progress=science_fit_progress,
                            coeff=science_fit_coeff,
                            linear=science_fit_linear,
                            weighted=science_fit_weighted,
                            filter_close=science_fit_filter_close,
                            display=science_fit_display,
                            savefig=science_fit_savefig,
                            filename=science_fit_filename,
                            spec_id=science_fit_spec_id,
                            stype='science')

            if science_refinefit1:
                onedspec[i].refine_fit(
                    polyfit_coeff=science_refinefit1_polyfit_coeff,
                    n_delta=science_refinefit1_n_delta,
                    refine=science_refinefit1_refine,
                    tolerance=science_refinefit1_tolerance,
                    method=science_refinefit1_method,
                    convergence=science_refinefit1_convergence,
                    robust_refit=science_refinefit1_robust_refit,
                    polydeg=science_refinefit1_polydeg,
                    display=science_refinefit1_display,
                    savefig=science_refinefit1_savefig,
                    filename=science_refinefit1_filename,
                    spec_id=science_refinefit1_spec_id,
                    stype='science')

            if science_refinefit2:
                onedspec[i].refine_fit(
                    polyfit_coeff=science_refinefit2_polyfit_coeff,
                    n_delta=science_refinefit2_n_delta,
                    refine=science_refinefit2_refine,
                    tolerance=science_refinefit2_tolerance,
                    method=science_refinefit2_method,
                    convergence=science_refinefit2_convergence,
                    robust_refit=science_refinefit2_robust_refit,
                    polydeg=science_refinefit2_polydeg,
                    display=science_refinefit2_display,
                    savefig=science_refinefit2_savefig,
                    filename=science_refinefit2_filename,
                    spec_id=science_refinefit2_spec_id,
                    stype='science')

            onedspec[i].apply_wavelength_calibration(stype='science')


#################################################################

##### #     #   # #   #
#     #     #   #  # #
##### #     #   #   #
#     #     #   #  # #
#     ##### ##### #   #

##### ##### #     ##### ####  ####  ##### ##### ##### ##### #   #
#     #   # #       #   #   # #   # #   #   #     #   #   # ##  #
#     ##### #       #   ####  ####  #####   #     #   #   # # # #
#     #   # #       #   #   # #  #  #   #   #     #   #   # #  ##
##### #   # ##### ##### ####  #   # #   #   #   ##### ##### #   #

#################################################################
try:
    sensitivity = np.load(sensitivity_filepath)
    for i in range(n_science):
        onedspec[i].add_sensitivity_itp(sensitivity)
    print('sensitivity_itp imported.')
    sensitivity_imported = True
except:
    if case == 2:
        print('sensitivity_itp cannot be imported or computed.')
    if case == 3:
        print('sensitivity_itp will be computed.')
    sensitivity_imported = False

if not sensitivity_imported:

    if case >= 3:
        standard_fits = fits.open(standard_filepath)[0]
        standard_flattened = standard_fits.data
        standard_header = standard_fits.header
        standard_northsouth = standard_header['NS']
        standard_binx = standard_header['BINX']

        if standard_pixel_list is None:
            standard_pixel_list = create_pixel_array(standard_northsouth,
                                                     standard_binx)

        if standard_atlas_relative_humidity is None:
            standard_atlas_relative_humidity = float(
                standard_header['HUMIDITY'])
        if standard_atlas_temperature is None:
            standard_atlas_temperature = float(
                standard_header['TAMBIENT']) + 273.15
        if standard_atlas_pressure is None:
            standard_atlas_pressure = float(standard_header['PRESSUR2'])

        if standard_atlas_user_relative_humidity is None:
            standard_atlas_user_relative_humidity = float(
                standard_header['HUMIDITY'])
        if standard_atlas_user_temperature is None:
            standard_atlas_user_temperature = float(
                standard_header['TAMBIENT']) + 273.15
        if standard_atlas_user_pressure is None:
            standard_atlas_user_pressure = float(standard_header['PRESSUR2'])

        # Make sure standard_noise, gain, seeing, exptime are None for the
        # above keywords to be used.
        standard_twodspec = spectral_reduction.TwoDSpec(
            standard_flattened,
            saxis=standard_saxis,
            spatial_mask=standard_spatial_mask,
            spec_mask=standard_spec_mask,
            flip=standard_flip,
            cosmicray=standard_cosmicray,
            cosmicray_sigma=standard_cosmicray_sigma,
            readnoise=standard_readnoise,
            gain=standard_gain,
            seeing=standard_seeing,
            exptime=standard_exptime,
            silence=standard_silence)

        standard_twodspec.ap_trace(
            nspec=standard_aptrace_nspec,
            nwindow=standard_aptrace_nwindow,
            spec_sep=standard_aptrace_spec_sep,
            resample_factor=standard_aptrace_resample_factor,
            rescale=standard_aptrace_rescale,
            scaling_min=standard_aptrace_scaling_min,
            scaling_max=standard_aptrace_scaling_max,
            scaling_step=standard_aptrace_scaling_step,
            percentile=standard_aptrace_percentile,
            tol=standard_aptrace_tol,
            polydeg=standard_aptrace_polydeg,
            ap_faint=standard_aptrace_ap_faint,
            display=standard_aptrace_display,
            renderer=standard_aptrace_renderer,
            return_jsonstring=standard_aptrace_jsonstring,
            save_iframe=standard_aptrace_iframe,
            open_iframe=standard_aptrace_open_iframe)

        standard_twodspec.ap_extract(
            apwidth=standard_apextract_apwidth,
            skysep=standard_apextract_skysep,
            skywidth=standard_apextract_skywidth,
            skydeg=standard_apextract_skydeg,
            optimal=standard_apextract_optimal,
            display=standard_apextract_display,
            renderer=standard_apextract_renderer,
            return_jsonstring=standard_apextract_jsonstring,
            save_iframe=standard_apextract_iframe,
            open_iframe=standard_apextract_open_iframe)

        onedspec[0].from_twodspec(standard_twodspec,
                                 pixel_list=standard_pixel_list,
                                 stype='standard')

    if case == 3:
        standard_arc_fits = fits.open(standard_arc_filepath)[0]
        standard_arc_flattened = standard_arc_fits.data
        standard_arc_header = standard_arc_fits.header
        standard_arc_northsouth = standard_arc_header['NS']
        standard_arc_binx = standard_arc_header['BINX']

        # Check arc binning matches other frames
        # GMOS-South didn't take any CuAr R400 in 1x1 binning since 2018
        if standard_arc_binx == standard_binx:
            onedspec[0].add_arc(np.flip(standard_arc_flattened),
                                stype='standard')
            onedspec[0].extract_arc_spec(
                display=standard_arc_spec_display,
                return_jsonstring=standard_arc_spec_jsonstring,
                renderer=standard_arc_spec_renderer,
                save_iframe=standard_arc_spec_iframe,
                open_iframe=standard_arc_spec_open_iframe,
                stype='standard')

        else:
            binx_ratio = standard_arc_binx / standard_binx
            spec = onedspec.standard_twodspec[0]

            len_trace = len(spec.trace)
            trace = np.nanmean(spec.trace) / binx_ratio
            trace_width = np.nanmean(spec.trace_sigma) * 3. / binx_ratio

            arc_trace = np.flip(standard_arc_flattened)[
                max(0, int(trace - trace_width -
                           1)):min(int(trace + trace_width), len_trace), :]
            arc_spec = np.nanmedian(arc_trace, axis=0)
            arc_spec_resampled = signal.resample(arc_spec, len_trace)

            onedspec[0].add_arc_spec(spec_id=0, arc_spec=arc_spec_resampled)

        onedspec[0].find_arc_lines(
            background=standard_findarc_background,
            percentile=standard_findarc_percentile,
            prominence=standard_findarc_prominence,
            distance=standard_findarc_distance,
            display=standard_findarc_display,
            refine=standard_findarc_refine,
            refine_window_width=standard_findarc_refine_window_width,
            return_jsonstring=standard_findarc_jsonstring,
            renderer=standard_findarc_renderer,
            save_iframe=standard_findarc_iframe,
            open_iframe=standard_findarc_open_iframe,
            stype='standard')
        onedspec[0].initialise_calibrator(
            pixel_list=standard_pixel_list,
            min_wavelength=standard_wavecal_min_wavelength,
            max_wavelength=standard_wavecal_max_wavelength,
            plotting_library=standard_wavecal_plotting_library,
            log_level=standard_wavecal_log_level,
            stype='standard')

        if standard_wavecal_polyfit:
            onedspec[0].add_polyfit(
                polyfit_coeff=standard_wavecal_polyfit_coeff,
                polyfit_type=standard_wavecal_polyfit_type,
                stype='standard')
        else:
            onedspec[0].set_fit_constraints(
                num_slopes=standard_constraints_num_slopes,
                range_tolerance=standard_constraints_range_tolerance,
                fit_tolerance=standard_constraints_fit_tolerance,
                polydeg=standard_constraints_polydeg,
                candidate_thresh=standard_constraints_candidate_thresh,
                linearity_thresh=standard_constraints_linearity_thresh,
                ransac_thresh=standard_constraints_ransac_thresh,
                num_candidates=standard_constraints_num_candidates,
                xbins=standard_constraints_xbins,
                ybins=standard_constraints_ybins,
                brute_force=standard_constraints_brute_force,
                polyfit_type=standard_constraints_polyfit_type,
                spec_id=standard_constraints_spec_id,
                stype='standard')

        if standard_atlas_user_supplied:
            if len(standard_atlas_user_elements) == 1:
                elements = [standard_atlas_user_elements
                            ] * len(standard_atlas_user_wavelengths)
            elif len(standard_atlas_user_elements) == len(
                    standard_atlas_user_wavelengths):
                elements = standard_atlas_user_elements
            else:
                raise ValueError('The science_atlas_user_elements should '
                                 'have length 1 or the same length as '
                                 'science_atlas_user_wavelengths')

            onedspec[0].load_user_atlas(
                elements=elements,
                wavelengths=standard_atlas_user_wavelengths,
                constrain_poly=standard_atlas_user_constrain_poly,
                vacuum=standard_atlas_user_vacuum,
                pressure=standard_atlas_user_pressure,
                temperature=standard_atlas_user_temperature,
                relative_humidity=standard_atlas_user_relative_humidity,
                stype='standard')
        else:
            onedspec[0].add_atlas(
                elements=standard_atlas_elements,
                min_atlas_wavelength=standard_atlas_min_atlas_wavelength,
                max_atlas_wavelength=standard_atlas_max_atlas_wavelength,
                min_intensity=standard_atlas_min_intensity,
                min_distance=standard_atlas_min_distance,
                vacuum=standard_atlas_vacuum,
                pressure=standard_atlas_pressure,
                temperature=standard_atlas_temperature,
                relative_humidity=standard_atlas_relative_humidity,
                constrain_poly=standard_atlas_constrain_poly,
                spec_id=standard_atlas_spec_id,
                stype='standard')

        onedspec[0].fit(sample_size=standard_fit_sample_size,
                        top_n=standard_fit_top_n,
                        max_tries=standard_fit_max_tries,
                        progress=standard_fit_progress,
                        coeff=standard_fit_coeff,
                        linear=standard_fit_linear,
                        weighted=standard_fit_weighted,
                        filter_close=standard_fit_filter_close,
                        display=standard_fit_display,
                        savefig=standard_fit_savefig,
                        filename=standard_fit_filename,
                        spec_id=standard_fit_spec_id,
                        stype='standard')

        if standard_refinefit1:
            onedspec[0].refine_fit(
                polyfit_coeff=standard_refinefit1_polyfit_coeff,
                n_delta=standard_refinefit1_n_delta,
                refine=standard_refinefit1_refine,
                tolerance=standard_refinefit1_tolerance,
                method=standard_refinefit1_method,
                convergence=standard_refinefit1_convergence,
                robust_refit=standard_refinefit1_robust_refit,
                polydeg=standard_refinefit1_polydeg,
                display=standard_refinefit1_display,
                savefig=standard_refinefit1_savefig,
                filename=standard_refinefit1_filename,
                spec_id=standard_refinefit1_spec_id,
                stype='standard')

        if standard_refinefit2:
            onedspec[0].refine_fit(
                polyfit_coeff=standard_refinefit2_polyfit_coeff,
                n_delta=standard_refinefit2_n_delta,
                refine=standard_refinefit2_refine,
                tolerance=standard_refinefit2_tolerance,
                method=standard_refinefit2_method,
                convergence=standard_refinefit2_convergence,
                robust_refit=standard_refinefit2_robust_refit,
                polydeg=standard_refinefit2_polydeg,
                display=standard_refinefit2_display,
                savefig=standard_refinefit2_savefig,
                filename=standard_refinefit2_filename,
                spec_id=standard_refinefit2_spec_id,
                stype='standard')

        onedspec[0].apply_wavelength_calibration(stype='standard')
        if n_science > 1:
            for i in range(1, n_science):
                onedspec[i].from_twodspec(standard_twodspec,
                                         pixel_list=standard_pixel_list,
                                         stype='standard')
                onedspec[i].wavecal_standard = onedspec[0].wavecal_standard
                onedspec[i].apply_wavelength_calibration(stype='standard')

# if the sensitivity_itp was not imported, it has to be computed
if case >= 3:

    if not sensitivity_imported:

        onedspec[0].load_standard(target=fluxcal_target,
                                  library=fluxcal_library,
                                  ftype=fluxcal_ftype,
                                  cutoff=fluxcal_cutoff)

        onedspec[0].compute_sensitivity(kind=sensecurve_kind,
                                        smooth=sensecurve_smooth,
                                        slength=sensecurve_slength,
                                        sorder=sensecurve_sorder,
                                        mask_range=sensecurve_mask_range,
                                        mask_fit_order=sensecurve_mask_fit_order,
                                        mask_fit_size=sensecurve_mask_fit_size)
        sensitivity_imported = True

        if sensitivity_filepath is not None:
            onedspec[0].save_sensitivity_itp(
                os.path.join(folder_path, sensitivity_filepath))

        if n_science > 1:
            for i in range(1, n_science):
                onedspec[i].fluxcal.spectrum_list_standard = onedspec[
                    0].fluxcal.spectrum_list_standard
                for j in range(science_aptrace_nspec):
                    onedspec[i].fluxcal.spectrum_list_science[
                        j].sensitivity_itp = onedspec[
                            0].fluxcal.spectrum_list_science[j].sensitivity_itp

for i in range(n_science):

    if case == 1:
        pass
    elif case == 3:
        onedspec[i].apply_flux_calibration(stype='science+standard')
    elif sensitivity_imported:
        if case == 3:
            onedspec[i].apply_flux_calibration(stype='science+standard')
        elif case == 2:
            onedspec[i].apply_flux_calibration(stype='science')
        else:
            pass
    else:
        warnings.warn('Flux calibration cannot be performed. Please check '
                      'the filepaths and the standard stars parameters.')

    # Save to FITS
    if "fits" in science_output_type:
        # Save as a FITS file
        onedspec[i].save_fits(output=fits_output,
                              filename=science_output_filepath[i],
                              stype='science',
                              overwrite=True)

    # Save to CSVs
    if "csv" in science_output_type:
        # save as CSV
        onedspec[i].save_csv(output=csv_output,
                             filename=science_output_filepath[i],
                             stype='science',
                             overwrite=True)

    if case >= 2:
        # Other save options:
        if save_png or save_png or save_svg or save_pdf or save_json or save_iframe:

            onedspec[i].inspect_reduced_spectrum(
                stype='science',
                wave_min=science_wavecal_min_wavelength,
                wave_max=science_wavecal_max_wavelength,
                renderer='default',
                width=width,
                height=height,
                filename=science_output_filepath[i],
                save_png=save_png,
                save_jpg=save_jpg,
                save_svg=save_svg,
                save_pdf=save_pdf,
                return_jsonstring=save_json,
                save_iframe=save_iframe,
                open_iframe=False)

if case >= 3:

    if "fits" in standard_output_type:
        # Save as a FITS file
        onedspec[0].save_fits(output=fits_output,
                              filename=standard_output_filepath,
                              stype='standard',
                              overwrite=True)

    if "csv" in standard_output_type:
        # save as CSV
        onedspec[0].save_csv(output=csv_output,
                             filename=standard_output_filepath,
                             stype='standard',
                             overwrite=True)

    if save_png_standard or save_png_standard or save_svg_standard or save_pdf_standard or save_json_standard or save_iframe_standard:

        onedspec[0].inspect_reduced_spectrum(
            stype='standard',
            wave_min=standard_wavecal_min_wavelength,
            wave_max=standard_wavecal_max_wavelength,
            renderer='default',
            width=width,
            height=height,
            filename=standard_output_filepath,
            save_png=save_png_standard,
            save_jpg=save_jpg_standard,
            save_svg=save_svg_standard,
            save_pdf=save_pdf_standard,
            return_jsonstring=save_json_standard,
            save_iframe=save_iframe_standard,
            open_iframe=False)
