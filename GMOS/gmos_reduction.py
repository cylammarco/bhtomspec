import os
import sys
import warnings
from copy import deepcopy
from time import time

from astropy.io import fits
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.lib import units
from scipy import ndimage
from svglib.svglib import svg2rlg
import yaml

from .gmos_fieldflattening import *

base_path = os.path.abspath(os.path.dirname(__file__))

# Get config file
try:
    params_path = os.path.join(base_path, sys.argv[1])
except:
    params_path = os.path.join(base_path, 'flattening_config.yaml')

if not os.path.isabs(params_path):
    params_path = os.path.abspath(params_path)

print('Reading parameters from ' + params_path + '.')

with open(params_path, 'r') as stream:
    params = yaml.safe_load(stream)

# Get the working directory and output directory
params_folder_path = os.path.dirname(params_path)
folder_path = os.path.join(params_folder_path, params['folder_path'])
output_folder_path = os.path.join(params_folder_path,
                                  params['output_folder_path'])

if not os.path.isabs(folder_path):
    folder_path = os.path.abspath(folder_path)

if output_folder_path is None:
    output_folder_path = folder_path

if not os.path.isabs(output_folder_path):
    output_folder_path = os.path.join(folder_path, output_folder_path)

if not os.path.exists(output_folder_path):
    os.mkdir(output_folder_path)

# Get all the parameters

# Hamamatsu detectors do not need dark subtraction
light_folder = params['light_folder']
flat_folder = params['flat_folder']
bias_folder = params['bias_folder']
arc_folder = params['arc_folder']
arc_flat_folder = params['arc_flat_folder']
arc_bias_folder = params['arc_bias_folder']

light_extension = params['light_extension']
flat_extension = params['flat_extension']
bias_extension = params['bias_extension']
arc_extension = params['arc_extension']
arc_flat_extension = params['arc_flat_extension']
arc_bias_extension = params['arc_bias_extension']

bias_master_filename = params['bias_master_filename']
arc_bias_master_filename = params['arc_bias_master_filename']

row_size = params['flat_row_size']
edge_size = params['flat_edge_size']
strip_size = params['flat_strip_size']

force_recreate_bias = params['force_recreate_bias']
save_bias = params['save_bias']
save_bias_format = params['save_bias_format']
overwrite_bias_image = params['overwrite_bias_image']

save_flattened_image = params['save_flattened_image']
save_flattened_image_format = params['save_flattened_image_format']
overwrite_flattened_image = params['overwrite_flattened_image']

create_fig = params['create_fig']
show_fig = params['show_fig']

# Does not work if create_fig is set to false
return_imgdata = params['diagnostic_pdf']

# Get the folder paths
if light_folder is None:
    light_folder = os.path.join(folder_path, 'light')

if flat_folder is None:
    flat_folder = os.path.join(folder_path, 'flat')

if bias_folder is None:
    bias_folder = os.path.join(folder_path, 'bias')

if arc_folder is None:
    arc_folder = os.path.join(folder_path, 'arc')

if arc_flat_folder is None:
    arc_flat_folder = os.path.join(folder_path, 'arc_flat')

if arc_bias_folder is None:
    arc_bias_folder = os.path.join(folder_path, 'arc_bias')

# If light or flat folders are empty, fail right away
light_path = os.listdir(light_folder)
if light_path == []:
    raise Error(light_folder + ' cannot be empty.')

flat_path = os.listdir(flat_folder)
if flat_path == []:
    raise Error(flat_folder + ' cannot be empty.')

'''
The process will be assignd as one of the nine scenarios.

* Bias on arc has basically no effect in the wavelength solution, we ignore
the bias if it is not specifically provided for the arcs.
* Flat on arc affects the quality of the find_peaks algorithm, so the flats
from the light frames will be used for crude reduction.

(1) Light + Flat

(2) Light + Flat + Bias

(3) Light + Flat + Arc (+ reuse Light Flat)

(4) Light + Flat + Bias + Arc (+ reuse Light Flat)

(5) Light + Flat + Arc + Arc Flat

(6) Light + Flat + Arc + Arc Bias (+ reuse Light Flat)

(7) Light + Flat + Arc + Arc Flat + Arc Bias

(8) Light + Flat + Bias + Arc + Arc Flat

(9) Light + Flat + Bias + Arc + Arc Bias (+ reuse Light Flat)

(10) Light + Flat + Bias + Arc + Arc Flat + Arc Bias

'''

# Preliminary assignment of the cases
# Folders are checked to exist, but contents are not guaranteed
case = 1

if os.path.exists(bias_folder):
    case = 2

if os.path.exists(arc_folder):
    case = 3

if os.path.exists(bias_folder) and os.path.exists(arc_folder):
    case = 4

if case == 3:
    if os.path.exists(arc_flat_folder):
        case = 5

    if os.path.exists(arc_bias_folder):
        case = 6

    if os.path.exists(arc_flat_folder) and os.path.exists(arc_bias_folder):
        case = 7

if case == 4:
    if os.path.exists(arc_flat_folder):
        case = 8

    if os.path.exists(arc_bias_folder):
        case = 9

    if os.path.exists(arc_flat_folder) and os.path.exists(arc_bias_folder):
        case = 10


# Get only the files with the right extension
light_path_new = []
for i, path_i in enumerate(light_path):
    if not os.path.isabs(path_i):
        path_i = os.path.abspath(os.path.join(light_folder, path_i))

    if os.path.splitext(path_i.lower())[-1][1:] in light_extension:
        if os.path.exists(path_i):
            light_path_new.append(path_i)
        else:
            raise ValueError(path_i + ' does not exist.')

light_path = light_path_new

# For all cases, flat frame must exist
flat_path_new = []
for i, path_i in enumerate(flat_path):
    if not os.path.isabs(path_i):
        path_i = os.path.abspath(os.path.join(flat_folder, path_i))

    if os.path.splitext(path_i.lower())[-1][1:] in flat_extension:
        if os.path.exists(path_i):
            flat_path_new.append(path_i)
        else:
            raise ValueError(path_i + ' does not exist.')

flat_path = flat_path_new


# If bias is available
if case in [2, 4, 8, 9, 10]:
    bias_northsouth = None
    # Load or create the bias master
    if os.path.isfile(
            os.path.join(output_folder_path, bias_master_filename +
                         '.npy')) and (not force_recreate_bias):
        print('I am loading the npy bias.')
        bias_master = np.load(
            os.path.join(output_folder_path, bias_master_filename + '.npy'))

        # assuming the bias_master is from the right telescope with the right
        # binning
        bias_binx = int(2048 / np.shape(bias_master)[2])
        bias_biny = int(4176 / np.shape(bias_master)[1])
        try:
            bias_master_filelist = np.loadtxt(
                os.path.join(output_folder_path,
                             bias_master_filename + '_filelist.txt'))
        except:
            bias_master_filelist = ['bias_master.npy']

    elif os.path.isfile(
            os.path.join(output_folder_path, bias_master_filename +
                         '.fits')) and (not force_recreate_bias):
        print('I am loading the FITS bias.')

        bias_master = fits.open(
            os.path.join(output_folder_path, bias_master_filename + '.fits'))[0]
        bias_master_filelist = []
        for i in bias_master[0].header:
            if i.startswith('BIAS'):
                bias_master_filelist.append(bias_master[0].header[i])

        # assuming the bias_master is from the right telescope with the right
        # binning
        bias_binx = int(2048 / np.shape(bias_master)[2])
        bias_biny = int(4176 / np.shape(bias_master)[1])
    else:
        try:
            print('I am creating the bias.')
            bias_fits, bias_binx, bias_biny, bias_northsouth = make_bias_master(
                bias_folder,
                bias_extension=bias_extension,
                save_bias=save_bias,
                save_bias_name=bias_master_filename,
                save_bias_format=save_bias_format,
                save_folder_path=output_folder_path,
                overwrite=overwrite_bias_image,
                create_fig=False,
                show_fig=False,
                return_imgdata=False)
            bias_master = bias_fits.data
            bias_header = bias_fits.header
            bias_master_filelist = []
            for i in bias_header:
                if i.startswith('BIAS'):
                    bias_master_filelist.append(bias_header[i])
        except:
            warnings.warn('Failed to creat the master bias. Reduction '
                'continues without bias correction.')
            if case == 2:
                case = 1
            if case == 4:
                case = 3
            if case == 8:
                case = 5
            if case == 9:
                case = 6
            if case == 10:
                case = 7

# Arc and thir flats and biases if different to the set for the light frames

# If arc is available
if case >=3:
    arc_path = os.listdir(arc_folder)
    for i, path_i in enumerate(arc_path):
        if not os.path.isabs(path_i):
            arc_path[i] = os.path.abspath(os.path.join(arc_folder, path_i))

    arc_path_new = []
    for i, path_i in enumerate(arc_path):
        if os.path.splitext(path_i.lower())[-1][1:] in arc_extension:
            if os.path.exists(path_i):
                arc_path_new.append(path_i)
            else:
                warnings.warn(path_i + ' does not exist.')
        else:
            pass

    arc_path = arc_path_new

    # Cases are reduced to 1 or 2 if arc frame is not avaialble
    if arc_path == []:
        warnings.warn(arc_folder + ' is either empty or does not contain '
            'files with the extensions: ' + arc_extension +  '.')

        if case in [3, 5, 6, 7]:
            case = 1
        if case == 4:
            case = 2
        if case >= 8:
            case = 2


# Check if arc flat is available
if case in [5, 7, 8, 10]:

    arc_flat_path = os.listdir(arc_flat_folder)

    arc_flat_path_new = []
    for i, path_i in enumerate(arc_flat_path):
        if not os.path.isabs(path_i):
            path_i = os.path.abspath(os.path.join(arc_flat_folder, path_i))
        if os.path.splitext(
                path_i.lower())[-1][1:] in arc_flat_extension:
            if os.path.exists(path_i):
                arc_flat_path_new.append(path_i)
            else:
                warnings.warn(path_i + ' does not exist.')
        else:
            pass

    arc_flat_path = arc_flat_path_new

    if arc_flat_path == []:
        warnings.warn(arc_flat_folder + ' is either empty or does not contain '
            'files with the extensions: ' + arc_flat_extension +  '. The '
            'arcs will be flattened with the science flat frames.')

        # If arc flats are not avaialbe, case is reduced
        if case == 5:
            case = 3
        if case == 7:
            case = 6
        if case == 8:
            case = 4
        if case == 10:
            case = 9


# If arc bias is available
if case in [6, 7, 9, 10]:

    arc_bias_path = os.listdir(arc_bias_folder)

    arc_bias_path_new = []
    for i, path_i in enumerate(arc_bias_path):
        if not os.path.isabs(path_i):
            path_i = os.path.abspath(os.path.join(arc_bias_folder, path_i))
        if os.path.splitext(
                path_i.lower())[-1][1:] in arc_bias_extension:
            arc_bias_path_new.append(path_i)
            if os.path.exists(path_i):
                arc_bias_path_new.append(path_i)
            else:
                warnings.warn(path_i + ' does not exist.')
        else:
            pass

    arc_bias_path = arc_bias_path_new

    if arc_bias_path == []:
        warnings.warn(arc_bias_folder + ' is either empty or does not contain '
            'files with the extensions: ' + arc_bias_extension +  '. The '
            'arcs will be flattened with the science flat frames.')

        # If arc bias are not avaialbe, case is reduced
        if case == 6:
            case = 3
        if case == 7:
            case = 5
        if case == 9:
            case = 4
        if case == 10:
            case = 8


# Reduction stars here

# Reconstruct the light frames
light_rc_frames = []
light_rc_bytesio_data = []
light_header = []
light_binx = []
light_biny = []
light_northsouth = []

for i, path_i in enumerate(light_path):
    light_temp = gmos_hamamatsu(fits.open(path_i),
                                create_fig=create_fig,
                                return_imgdata=return_imgdata)
    light_rc_frames.append(light_temp[0])
    light_binx.append(light_temp[2])
    light_biny.append(light_temp[3])
    light_northsouth.append(light_temp[4])
    light_header.append(light_temp[5])
    # Store the diagnostic plots if return_imgdata is True
    if return_imgdata:
        light_rc_bytesio_data.append(light_temp[-1])

if len(set(light_northsouth)) > 1:
    raise ValueError(
        'The light frames contain a mix of GMOST North and South.')
else:
    light_northsouth = light_northsouth[0]

if len(set(light_binx)) > 1:
    # larger binning -> lower resolution
    binx_min = min(light_binx)
    resample_idx = np.where(light_binx > binx_min)[0]
    for i in resample_idx:
        warnings.warn('Light frames have different binnings, the lower '
                      'resolution images are UPsampled.')
        light_rc_frames[i] = ndimage.zoom(light_rc_frames[i],
                                          light_binx[i] / binx_min)
    light_binx = binx_min
    light_biny = binx_min
else:
    light_binx = light_binx[0]
    light_biny = light_biny[0]


if case in [2, 4, 8, 9, 10]:
    # If the binning of the bias frame does not match that of the light frames,
    # up/down sample it and throw a warning
    if bias_binx != light_binx:
        warnings.warn(
            'Light frame and bias frame have different binning, '
            'the bias is resampled to match the light frame. Total count may '
            'not be conserved.')
        bias_master = ndimage.zoom(bias_master, bias_binx / light_binx)

    # If the bias is not known to be from the GMOS north or south, assume it's
    # the right one
    if bias_northsouth is None:
        bias_northsouth = light_northsouth

    bias_pixels = create_pixel_array(bias_northsouth, bias_binx)

    bias_response_bytesio_data = plot_bias_response(bias_master,
                                                    pixels=bias_pixels,
                                                    show_fig=show_fig,
                                                    return_imgdata=return_imgdata)

# Reconstruct the flat frames
flat_rc_frames = []
flat_rc_bytesio_data = []

for i, path_i in enumerate(flat_path):
    flat_temp = gmos_hamamatsu(fits.open(path_i),
                               create_fig=create_fig,
                               return_imgdata=return_imgdata)
    flat_rc_frames.append(flat_temp[0])
    flat_binx = flat_temp[2]
    flat_biny = flat_temp[3]
    flat_northsouth = flat_temp[4]

    if flat_binx != light_binx:
        warnings.warn(
            'Light frame and flat frame have different binning, '
            'the flat is resampled to match the light frame. Total count may '
            'not be conserved.')
        flat_rc_frames[i] = ndimage.zoom(flat_rc_frames[i],
                                         flat_binx / light_binx)

    # Store the diagnostic plots if return_imgdata is True
    if return_imgdata:
        flat_rc_bytesio_data.append(flat_temp[-1])


# Reconstruct the arc frames
arc_rc_frames = []
arc_rc_bytesio_data = []
arc_header = []
arc_binx = []
arc_biny = []
if case >= 3:

    for i, path_i in enumerate(arc_path):
        arc_temp = gmos_hamamatsu(fits.open(path_i),
                                  create_fig=create_fig,
                                  return_imgdata=return_imgdata)
        arc_rc_frames.append(arc_temp[0])
        arc_binx.append(arc_temp[2])
        arc_biny.append(arc_temp[3])
        arc_northsouth = arc_temp[4]
        arc_header.append(arc_temp[5])
        # Store the diagnostic plots if return_imgdata is True
        if return_imgdata:
            arc_rc_bytesio_data.append(arc_temp[-1])

    if len(set(arc_binx)) > 1:
        # larger binning -> lower resolution
        binx_min = min(arc_binx)
        resample_idx = np.where(arc_binx > binx_min)[0]
        arc_rc_frames_temp = []
        for i in resample_idx:
            warnings.warn('Arc frames have different binnings, the lower '
                          'resolution images are UPsampled.')
            arc_rc_frames_temp.append(
                ndimage.zoom(arc_rc_frames[i], arc_binx[i] / binx_min))
        arc_rc_frames = arc_rc_frames_temp
        arc_binx = binx_min
        arc_biny = binx_min
    else:
        arc_binx = arc_binx[0]
        arc_biny = arc_biny[0]


# Reconstruct the arc flat frames
arc_flat_rc_frames = []
arc_flat_rc_bytesio_data = []
if case in [5, 7, 8, 10]:

    for i, path_i in enumerate(arc_flat_path):
        arc_flat_temp = gmos_hamamatsu(fits.open(path_i),
                                       create_fig=create_fig,
                                       return_imgdata=return_imgdata)
        arc_flat_rc_frames.append(flat_temp[0])
        arc_flat_binx = flat_temp[2]
        arc_flat_biny = flat_temp[3]
        arc_flat_northsouth = flat_temp[4]

        # Store the diagnostic plots if return_imgdata is True
        if return_imgdata:
            arc_flat_rc_bytesio_data.append(arc_flat_temp[-1])

    if arc_flat_binx != arc_binx:
        arc_flat_rc_frames = []
        for i in range(len(flat_rc_frames)):
            warnings.warn(
                'Light frame and flat frame have different binning, '
                'the flat is resampled to match the light frame. Total count may '
                'not be conserved.')
            arc_flat_rc_frames.append(
                ndimage.zoom(flat_rc_frames[i], arc_flat_binx / arc_binx))
        arc_flat_rc_frames = np.array(arc_flat_rc_frames)
    else:
        arc_flat_rc_frames = flat_rc_frames

if case in [6, 7, 9, 10]:
    # Load or create the arc bias master
    if os.path.isfile(
            os.path.join(output_folder_path, arc_bias_master_filename +
                         '.npy')) and (not force_recreate_bias):
        print('I am loading the npy arc bias.')
        arc_bias_master = np.load(
            os.path.join(output_folder_path,
                         arc_bias_master_filename + '.npy'))

        # assuming the arc_bias_master is from the right telescope with the right
        # binning
        arc_bias_northsouth = light_northsouth
        arc_bias_binx = int(2048 / np.shape(arc_bias_master)[2])
        arc_bias_biny = int(4176 / np.shape(arc_bias_master)[1])
        arc_bias_master_imported = True
        try:
            arc_bias_master_filelist = np.loadtxt(
                os.path.join(output_folder_path,
                             arc_bias_master_filename + '_filelist.txt'))
        except:
            arc_bias_master_filelist = ['arc_bias_master.npy']

    elif os.path.isfile(
            os.path.join(output_folder_path, arc_bias_master_filename +
                         '.fits')) and (not force_recreate_bias):
        print('I am loading the FITS arc bias.')

        arc_bias_master_fits = fits.open(
            os.path.join(output_folder_path,
                         arc_bias_master_filename + '.fits'))[0]
        arc_bias_master_filelist = []
        for i in arc_bias_master[0].header:
            if i.startswith('BIAS'):
                arc_bias_master_filelist.append(
                    arc_bias_master_fits[0].header[i])

        arc_bias_master = arc_bias_master_fits.data

        # assuming the bias_master is from the right telescope with the right
        # binning
        arc_bias_northsouth = light_northsouth
        arc_bias_binx = int(2048 / np.shape(arc_bias_master)[2])
        arc_bias_biny = int(4176 / np.shape(arc_bias_master)[1])
        arc_bias_master_imported = True
    else:
        print('I am creating the arc bias.')
        arc_bias_fits, arc_bias_binx, arc_bias_biny, arc_bias_northsouth = make_bias_master(
            arc_bias_folder,
            bias_extension=arc_bias_extension,
            save_bias=save_bias,
            save_bias_name=arc_bias_master_filename,
            save_bias_format=save_bias_format,
            save_folder_path=output_folder_path,
            overwrite=overwrite_bias_image,
            create_fig=False,
            show_fig=False,
            return_imgdata=False)
        arc_bias_master = arc_bias_fits.data
        arc_bias_header = arc_bias_fits.header
        arc_bias_master_filelist = []
        for i in arc_bias_header:
            if i.startswith('BIAS'):
                arc_bias_master_filelist.append(arc_bias_header[i])
        arc_bias_master_imported = True

    if arc_bias_binx != arc_binx:
        warnings.warn(
            'arc_bias frame and arc frame have different binning, the '
            'arc_bias is resampled to match the arc frame. Total count may '
            'not be conserved.')
        arc_bias_master = ndimage.zoom(arc_bias_master, arc_bias_binx / arc_binx)

    arc_bias_pixels = create_pixel_array(arc_bias_northsouth, arc_bias_binx)

    arc_bias_response_bytesio_data = plot_bias_response(
        arc_bias_master,
        pixels=arc_bias_pixels,
        show_fig=show_fig,
        return_imgdata=return_imgdata)


# Bias subtraction
if case in [2, 4, 8, 9, 10]:
    for frame_i in light_rc_frames:
        frame_i -= bias_master

    for frame_i in flat_rc_frames:
        frame_i -= bias_master


# Mean combine the flat frames
# Get the relative detector response from the flats, correcting for both
# CCD response and vignetting across the focal plane
flat = np.nanmean(flat_rc_frames, axis=0)
flat_normed, flat_normalisation, flat_stacked_bytesio_data = normalise_flat(
    flat=flat,
    binx=flat_binx,
    biny=flat_biny,
    northsouth=flat_northsouth,
    row_size=row_size,
    edge_size=edge_size,
    strip_size=strip_size,
    create_fig=create_fig,
    show_fig=show_fig,
    return_imgdata=return_imgdata)

# Flatten the light frames
light_flattened = []
light_flattened_bytesio_data = []
for frame_i in light_rc_frames:
    light_temp = flatten_image(image=frame_i,
                               flat=flat_normed,
                               binx=flat_binx,
                               biny=flat_biny,
                               normalisation=flat_normalisation,
                               create_fig=create_fig,
                               show_fig=show_fig,
                               return_imgdata=return_imgdata)
    if return_imgdata:
        light_flattened.append(light_temp[0])
        light_flattened_bytesio_data.append(light_temp[1])
    else:
        light_flattened.append(light_temp)


# the arc_flat is also referenced to the flat if it is not avaialble, so the
# else case was already bias-subtracted above.
if case in [5, 7, 8, 10]:

    if case in [8, 10]:
        for frame_i in arc_flat_rc_frames:
            frame_i -= arc_bias_master

    # Mean combine the flat frames
    # Get the relative detector response from the flats, correcting for both
    # CCD response and vignetting across the focal plane
    arc_flat = np.nanmean(arc_flat_rc_frames, axis=0)
    arc_flat_normed, arc_flat_normalisation, arc_flat_stacked_bytesio_data = normalise_flat(
        flat=arc_flat,
        binx=arc_flat_binx,
        biny=flat_biny,
        northsouth=arc_flat_northsouth,
        row_size=row_size,
        edge_size=edge_size,
        strip_size=strip_size,
        create_fig=create_fig,
        show_fig=show_fig,
        return_imgdata=return_imgdata)
# if the arc_flat is not available, the flat is used instead.
elif case in [3, 4, 6, 9]:
    arc_flat_normed = []
    arc_flat_normalisation = flat_normalisation
    arc_flat_stacked_bytesio_data = flat_stacked_bytesio_data

    # Make a deep copy if it is a reference to flat_rc_frames
    for i in range(len(flat_normed)):
        warnings.warn(
            'Arc frame and arc flat frame have different binning, '
            'the flat is resampled to match the light frame. Total count may '
            'not be conserved.')
        # use flat_binx instead of arc_flat_binx
        arc_flat_binx = flat_binx
        arc_flat_biny = flat_biny
        arc_flat_normed.append(
            ndimage.zoom(flat_normed[i], arc_flat_binx / arc_binx))
    arc_flat_normed = np.array(arc_flat_normed)
else:
    pass


# Flatten the arc frames
arc_flattened = []
arc_flattened_bytesio_data = []

if case >= 3:

    for frame_i in arc_rc_frames:
        arc_temp = flatten_image(image=frame_i,
                                 flat=arc_flat_normed,
                                 binx=arc_flat_binx,
                                 biny=arc_flat_biny,
                                 normalisation=arc_flat_normalisation,
                                 create_fig=create_fig,
                                 show_fig=show_fig,
                                 return_imgdata=return_imgdata)
        if return_imgdata:
            arc_flattened.append(arc_temp[0])
            arc_flattened_bytesio_data.append(arc_temp[1])
        else:
            arc_flattened.append(arc_temp)

    frame_path = light_path + arc_path
    frame_header = light_header + arc_header
    frame_flattened = light_flattened + arc_flattened
else:

    frame_path = light_path
    frame_header = light_header
    frame_flattened = light_flattened

# Save as npy or FITS
if save_flattened_image:

    for i, (frame_i, name_i) in enumerate(zip(frame_flattened, frame_path)):

        flattened_image = fits.PrimaryHDU(frame_i, header=fits.Header())

        # Add the names of the frames to header
        #flattened_image.header['COMMENT'] = "The frames."
        for i in range(len(flat_path)):
            flattened_image.header.set(
                keyword='FLAT' + str(i + 1),
                value=flat_path[i].split('/')[-1].split('.')[0],
                comment='Flat frame ' + str(i + 1))
        # Add other keywords
        if (i < len(light_path)) and (case in [2, 4, 8, 9, 10]):
            for i in range(len(bias_master_filelist)):
                for i, file_i in enumerate(bias_master_filelist):
                    flattened_image.header.set(keyword='bias' + str(i + 1),
                                               value=file_i,
                                               comment='Bias frame ' + str(i + 1))
        if (i >= len(light_path)) and (case in [6, 9, 10]):
            for i in range(len(arc_bias_master_filelist)):
                for i, file_i in enumerate(arc_bias_master_filelist):
                    flattened_image.header.set(keyword='bias' + str(i + 1),
                                               value=file_i,
                                               comment='Bias frame ' + str(i + 1))
        flattened_image.header.set(
            keyword='SCLIP',
            value='3',
            comment='Number of sigma used for outlier clipping.')
        flattened_image.header.set(keyword='NS',
                                   value=light_northsouth,
                                   comment='GMOS North or South.')
        flattened_image.header.set(
            keyword='BINX',
            value=light_binx,
            comment='Binning in the spectral direction.')
        flattened_image.header.set(keyword='BINY',
                                   value=light_biny,
                                   comment='Binning in the spatial direction.')
        flattened_image.header.set(keyword='PRESSUR2',
                                   value=frame_header[i]['PRESSUR2'],
                                   comment='Pressure (Pa).')
        flattened_image.header.set(keyword='TAMBIENT',
                                   value=frame_header[i]['TAMBIENT'],
                                   comment='Ambient temperature (C).')
        flattened_image.header.set(keyword='HUMIDITY',
                                   value=frame_header[i]['HUMIDITY'],
                                   comment='Relative humidity (0-100%).')

        flattened_image.writeto(os.path.join(
            output_folder_path,
            name_i.split('/')[-1].split('.')[0] + '_flattened.fits'),
                                overwrite=overwrite_flattened_image)

# Make PDF, one for each light frame
if return_imgdata:
    frame_list = zip(light_rc_bytesio_data + arc_rc_bytesio_data,
                     light_flattened_bytesio_data + arc_flattened_bytesio_data)
    filenames = frame_path

    flat_drawing = svg2rlg(flat_stacked_bytesio_data)

    for i, frame_i in enumerate(frame_list):

        if (i < len(light_path)) and (case in [2, 4, 8, 9, 10]):
            bias_drawing = svg2rlg(bias_response_bytesio_data)

        if (i >= len(light_path)) and (case in [6, 9, 10]):
            arc_bias_drawing = svg2rlg(arc_bias_response_bytesio_data)

        rc_drawing = svg2rlg(frame_i[0])
        flattened_drawing = svg2rlg(frame_i[1])

        outfile_name = filenames[i].split('/')[-1].split(
            '.')[0] + '_diagnostic_plot.pdf'
        outfile_path = os.path.join(output_folder_path, outfile_name)

        # If the file already exists, output filename with timestamp
        if os.path.isfile(outfile_path):
            tmp = outfile_name.split('.')
            tmp[-2] = tmp[-2] + "_" + str(time()).split('.')[0]
            outfile_name = ".".join(tmp)

        c = canvas.Canvas(outfile_path)
        c.scale(0.625, 0.625)
        c.setFont("Helvetica", 24)
        c.drawString(150 * units.mm, 460 * units.mm, outfile_name)
        renderPDF.draw(rc_drawing, c, 10 * units.mm, 300 * units.mm)
        if (i < len(light_path)) and (case in [2, 4, 8, 9, 10]):
            renderPDF.draw(bias_drawing, c, 10 * units.mm, 220 * units.mm)
        if (i >= len(light_path)) and (case in [6, 9, 10]):
            renderPDF.draw(arc_bias_drawing, c, 10 * units.mm, 220 * units.mm)
        renderPDF.draw(flat_drawing, c, 10 * units.mm, 150 * units.mm)
        renderPDF.draw(flattened_drawing, c, 10 * units.mm, 10 * units.mm)
        c.showPage()
        c.save()
