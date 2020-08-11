import os
from io import BytesIO

from astropy.io import fits
from astropy.stats import sigma_clip
from matplotlib import pyplot as plt
from scipy import optimize
import numpy as np

# All pixel values here are unbinned
# b for binning
# n for north
# s for south
# rc for reconstructed
# ro for readout
# os for overscan
hamamatsu_width = 2048
hamamatsu_height = 4176
hamamatsu_exposed_width = 512
hamamatsu_unused_height = 48
hamamatsu_os_width = 32

# chip gaps
hamamatsu_n_gap = 67
hamamatsu_s_gap = 61

# min and max for the central unobscured part of the CCDs
hamamatsu_middle_chip_min = 1450
hamamatsu_middle_chip_max = 2750


def create_pixel_array(northsouth, binning):
    '''
    The GMOS longslit spectrum spreads over three CCDs. This function creates
    the pixel list corresponding the binned pixels adjusted for the chip gaps,
    which can lead to a non-integer shift in the number of pixels.

    Parameters
    ----------
    northsout: str
        Indicate whether the Gemini north or south, because the chip gaps are
        different.
    binning : numeric
        The binning factor

    Returns
    -------
    pixels: numpy.array
        The pixel values adjusted for chip gaps for the corresponding pixel.

    '''

    binned_width = hamamatsu_width / binning

    if northsouth == 'north':
        gap_binned_width = hamamatsu_n_gap / binning
    elif northsouth == 'south':
        gap_binned_width = hamamatsu_s_gap / binning
    else:
        raise ValueError('Please choose from "north" or "south".')

    pixels = np.concatenate(
        (np.arange(binned_width),
         np.arange(binned_width + gap_binned_width,
                   binned_width * 2 + gap_binned_width),
         np.arange(binned_width * 2 + gap_binned_width * 2,
                   binned_width * 3 + gap_binned_width * 2)))

    return pixels


def gmos_hamamatsu(fits_data,
                   create_fig=True,
                   show_fig=False,
                   show_gap=True,
                   show_os=True,
                   return_imgdata=False):
    '''
    This function reconstruct the multiple extension FITS file into a flat
    image. The chip gaps are not relavent here. The overscan signals are
    subtracted during the image reconstruction.

    Note: ROI for CENTRAL SPECTRUM and CCD2 maps to a full frame.

    Parameters
    ----------
    fits_data: astropy.io.fits
        The multiple extension data file. The primary HDU does not hold data.
    binning: numeric
        The binning factor, which should be an integer, but data type is not
        important.
    create_fig: boolean
        Set to True to create a figure with matplotlib in the backend.
    show_fig: boolean
        Set to True to view the figure.
    show_gap: boolean
        Set to True to display (arbitrary) spacing between the chips.
    show_os: boolean
        Set to True to display the overscan region, which is marked by the
        black dashed box.
    return_imgdata: boolean
        Set to True to return an svg image as a BytesIO object.

    Returns
    -------
    rc_image: numpy.array
        The 3 reconstructed images for each of the CCDs.
    os: numpy.array
        The collapsed overscan 1D numpy array for each of the 12 readout
        regions.
    imgdata: BytesIO
        BytesIO object holding the svg image.

    '''

    # Get binning
    binx, biny = np.array(fits_data[1].header['CCDSUM'].split()).astype('int')
    northsouth = fits_data[0].header['TELESCOP'].split('-')[1].lower()

    # Check if central spectrum only:
    roi_y_size = int(fits_data[0].header['DETRO1YS'])  # iraf starts from 1
    roi_y_start = int(fits_data[0].header['DETRO1Y'])

    # Note that this ignores the chip gaps
    rc_width = int(hamamatsu_width / binx)
    rc_height = int(hamamatsu_height / biny)
    exposed_width_binned = int(hamamatsu_exposed_width / binx)
    unused_row = int(hamamatsu_unused_height / biny)

    # empty 2D array to receive the 3 reconstructed frames
    rc_image1 = np.zeros((rc_height, rc_width))
    rc_image2 = np.zeros((rc_height, rc_width))
    rc_image3 = np.zeros((rc_height, rc_width))

    # overscan
    os = []
    pfit = []

    if roi_y_size == 1024 / biny:
        roi_y_start = int(roi_y_start / biny) - 1
    else:
        roi_y_start -= 1

    if roi_y_size == 4224 / biny:
        roi_y_size -= unused_row
        start_row = unused_row
    else:
        start_row = 0

    roi_y_end = roi_y_start + roi_y_size
    end_row = start_row + roi_y_size

    for i in range(1, 13):
        start = int((i - 1) % 4 * exposed_width_binned)
        end = int(start + exposed_width_binned)
        data = fits_data[i].data
        header = fits_data[i].header

        if (i < 5):
            image = rc_image1
        elif (i >= 5) and (i < 9):
            image = rc_image2
        else:
            image = rc_image3

        # odd frames have overscan on the right, even frames have it on the left
        if i % 2 == 1:
            image[roi_y_start:roi_y_end,
                  start:end] = data[start_row:end_row, :exposed_width_binned]
            # correct with overscan
            os.append(np.mean(data[start_row:, exposed_width_binned:], axis=1))
            pfit.append(
                np.polyfit(np.arange(len(os[i - 1]))[20:-20],
                           os[i - 1][20:-20],
                           deg=15))
        else:
            image[roi_y_start:roi_y_end,
                  start:end] = data[start_row:end_row, -exposed_width_binned:]
            # correct with overscan
            os.append(np.mean(data[start_row:, :-exposed_width_binned],
                              axis=1))
            pfit.append(
                np.polyfit(np.arange(len(os[i - 1]))[20:-20],
                           os[i - 1][20:-20],
                           deg=15))

        image[roi_y_start:roi_y_end,
              start:end] -= np.polyval(pfit[i - 1],
                                       np.arange(len(os[i - 1])))[:, None]

    if not create_fig:
        return (rc_image1, rc_image2, rc_image3
                ), np.array(os), binx, biny, northsouth, fits_data[0].header

    else:
        if show_gap:
            gridspec = dict(
                wspace=0.0,
                left=0.07,
                right=0.99,
                top=0.99,
                bottom=0.05,
                width_ratios=[1, 1, 1, 1, 0.2, 1, 1, 1, 1, 0.2, 1, 1, 1, 1])
            fig, axes = plt.subplots(1,
                                     14,
                                     figsize=(10, 4.5),
                                     gridspec_kw=gridspec,
                                     sharex=True,
                                     sharey=True)
        else:
            fig, axes = plt.subplots(1,
                                     12,
                                     figsize=(10, 4.5),
                                     sharex=True,
                                     sharey=True)
        plt.gcf().set_size_inches(10, 4.5)

        # i for FITS index
        # j for axes index adjusted for empty axes to show the gaps
        for i in range(1, 13):
            j = i - 1

            if show_gap:
                if i >= 5 and i <= 8:
                    j += 1
                if i >= 9:
                    j += 2

            if i % 2 == 1:
                os_start = exposed_width_binned
                os_end = exposed_width_binned + hamamatsu_os_width
                im_start = 0
                im_end = exposed_width_binned
            else:
                os_start = 0
                os_end = hamamatsu_os_width
                im_start = hamamatsu_os_width
                im_end = hamamatsu_os_width + exposed_width_binned

            data = fits_data[i].data

            if show_os:
                # The 0.5 is for clerity purpose
                os_corners_x = [
                    os_start + 0.5, os_end - 0.5, os_end - 0.5, os_start + 0.5,
                    os_start + 0.5
                ]
                os_corners_y = [
                    0.5, 0.5, rc_height + start_row - 0.5,
                    rc_height + start_row - 0.5, 0.5
                ]
                axes[j].imshow(np.log10(data - np.nanpercentile(data, 5)),
                               aspect='auto',
                               origin='lower')
                axes[j].plot(os_corners_x,
                             os_corners_y,
                             color='black',
                             ls="--",
                             lw=1.5)
                axes[j].set_xlim(0., exposed_width_binned + hamamatsu_os_width)
                axes[j].set_ylim(0., roi_y_size)
            else:
                axes[j].imshow(
                    np.log10(data[:, im_start:im_end] -
                             np.nanpercentile(data[:, im_start:im_end], 5)),
                    aspect='auto',
                    origin='lower')

        # Show aribitrary chip gaps
        if show_gap:
            axes[4].set_visible(False)
            axes[9].set_visible(False)
            axes[6].set_title('Raw Image')
            axes[6].set_xlabel('Spectral Direction / Pixel')
            axes[0].set_ylabel('Spatial Direction / Pixel')
        else:
            axes[5].set_title('Raw Image')
            axes[5].set_xlabel('Spectral Direction / Pixel')
            axes[0].set_ylabel('Spatial Direction / Pixel')
            plt.subplots_adjust(wspace=0)
            fig.tight_layout()

        if show_fig:
            plt.show()

        if return_imgdata:
            imgdata = BytesIO()
            fig.savefig(imgdata, format='svg')
            imgdata.seek(0)  # rewind the data
            plt.close()
            return (rc_image1, rc_image2, rc_image3), np.array(
                os), binx, biny, northsouth, fits_data[0].header, imgdata
        else:
            plt.close()
            return (rc_image1, rc_image2, rc_image3), np.array(
                os), binx, biny, northsouth, fits_data[0].header


def make_bias_master(folder_path,
                     bias_extension=['fits'],
                     save_bias=True,
                     save_bias_name='bias_master',
                     save_bias_format=['npy'],
                     save_folder_path=None,
                     overwrite=True,
                     create_fig=True,
                     show_fig=False,
                     return_imgdata=False):
    '''
    Meancombine the bias frames.

    Parameters
    ----------
    folder_path: str
        The folder path containing the bias data only.
    save_bias: boolean
        Set to True to save the bias as an numpy array
    save_folder_path: str
        Destination to save the bias_master. Default is to save at the same
        folder holding the individual bias frames
    create_fig: boolean
        Set to True to create a figure with matplotlib in the backend.
    show_fig: boolean
        Set to True to view the figure.
    return_imgdata: boolean
        Set to True to return an svg image as a BytesIO object.

    Returns
    -------
    bias_master: numpy.array
        The combined bias frame.
    imgdata: BytesIO
        BytesIO object holding the svg image.

    '''

    bias_list_temp = os.listdir(folder_path)
    bias_list = []
    bias_rc = []

    for bias_path in bias_list_temp:
        if os.path.splitext(bias_path.lower())[-1][1:] in bias_extension:
            bias_list.append(bias_path)
            temp = fits.open(os.path.join(folder_path, bias_path))
            bias_frame, _, binx, biny, northsouth, bias_header = gmos_hamamatsu(
                temp, create_fig=False)
            bias_rc.append(bias_frame)

    bias_master = np.nanmean(sigma_clip(bias_rc,
                                        cenfunc=np.nanmean,
                                        axis=0,
                                        masked=False),
                             axis=0)

    # Put the reduced data in FITS format with an image header
    bias_fits = fits.PrimaryHDU(bias_master, header=fits.Header())

    # Add the names of all the biad frames to header
    bias_fits.header['COMMENT'] = "Bias frames used are listed below."
    if len(bias_list) > 0:
        for i in range(len(bias_list)):
            bias_fits.header.set(keyword='bias' + str(i + 1),
                                 value=bias_list[i],
                                 comment='Bias frame ' + str(i + 1))
    # Add other keywords
    bias_fits.header.set(keyword='COMBTYPE',
                         value='NANMEAN',
                         comment='Type of image combine.')
    bias_fits.header.set(keyword='SCLIP',
                         value='3',
                         comment='Number of sigma used for outlier clipping.')

    del bias_rc

    if save_bias:

        # Save as npy
        if 'npy' in save_bias_format:
            if save_folder_path is None:
                np.save(os.path.join(folder_path, save_bias_name + '.npy'),
                        bias_master)
                np.savetxt(s.path.join(folder_path,
                                       save_bias_name + '_filelist.txt'),
                           bias_list,
                           fmt="%s")
            else:
                np.save(
                    os.path.join(save_folder_path, save_bias_name + '.npy'),
                    bias_master)
                np.savetxt(os.path.join(save_folder_path,
                                        save_bias_name + '_filelist.txt'),
                           bias_list,
                           fmt="%s")

        # Save as FITS files
        if ('fits' or 'fit') in save_bias_format:

            if save_folder_path is None:

                bias_fits = fits.PrimaryHDU(biasfits)
                # Save file to disk
                bias_fits.writeto(os.path.join(folder_path,
                                               save_bias_name + '.fits'),
                                  overwrite=overwrite)

            else:
                bias_fits.writeto(os.path.join(save_folder_path,
                                               save_bias_name + '.fits'),
                                  overwrite=overwrite)

    if create_fig:
        pixels = create_pixel_array(northsouth, binx)
        if return_imgdata:
            imgdata = plot_bias_response(bias_master,
                                         pixels=pixels,
                                         show_fig=show_fig,
                                         return_imgdata=return_imgdata)
        else:
            plot_bias_response(bias_master,
                               pixels=pixels,
                               show_fig=show_fig,
                               return_imgdata=return_imgdata)

    return bias_fits, binx, biny, northsouth


def plot_bias_response(bias_master,
                       pixels=None,
                       show_fig=False,
                       return_imgdata=False):
    '''
    Plot the average bias signal at the central rows.

    Parameters
    ----------
    bias_master: numpy.array
        The combined bias frame.
    pixels: numpy.array
        The pixel values adjusted for chip gaps for the corresponding pixel.
    show_fig: boolean
        Set to True to view the figure.
    return_imgdata: boolean
        Set to True to return an svg image as a BytesIO object.

    Returns
    -------
    imgdata: BytesIO
        BytesIO object holding the svg image.
    '''

    if pixels is None:
        pixels = np.arange(len(bias_master[0][0]) * 3)
    mid_row = int(len(bias_master[0]) / 2)
    width = 100
    fig = plt.figure()
    plt.gcf().set_size_inches(10, 2.3)
    plt.plot(pixels[:int(len(pixels) / 3)],
             np.median(bias_master[0][mid_row - width:mid_row + width],
                       axis=0),
             color='red',
             label='Red CCD')
    plt.plot(pixels[int(len(pixels) / 3):int(len(pixels) / 3 * 2)],
             np.median(bias_master[1][mid_row - width:mid_row + width],
                       axis=0),
             color='green',
             label='Green CCD')
    plt.plot(pixels[int(len(pixels) / 3 * 2):],
             np.median(bias_master[2][mid_row - width:mid_row + width],
                       axis=0),
             color='blue',
             label='Blue CCD')
    plt.xlim(min(pixels), max(pixels))
    plt.ylim(np.nanpercentile(bias_master, 15),
             np.nanpercentile(bias_master, 95))
    plt.title('Bias Response')
    plt.xlabel('Spectral Direction / Pixel')
    plt.ylabel('ADU')
    plt.grid()
    plt.tight_layout()

    if show_fig:
        plt.show()

    if return_imgdata:
        imgdata = BytesIO()
        fig.savefig(imgdata, format='svg')
        imgdata.seek(0)  # rewind the data
        plt.close()
        return imgdata
    else:
        plt.close()


def linear_fit(guess, y1, y2, gap):
    '''
    Minimisation function to align the two-part linear function to a single
    straight line.
    '''

    m, norm, c = guess
    ratio = np.mean(y2) / np.mean(y1)
    if (norm > ratio * 2.) or (norm < ratio * 0.5) or (c < y1[0] * 0.5) or (
            c > y1[0] * 2.0):
        return np.inf
    x1 = np.arange(len(y1))
    x2 = np.arange(len(y2)) + len(y1) + gap
    model1 = m * x1 + c
    model2 = (m * x2 + c) * norm
    diff1 = model1 - y1
    diff2 = model2 - y2
    diff = np.concatenate((diff1, diff2))
    return np.sum(diff**2.)


def normalise_flat(flat,
                   binx,
                   biny,
                   northsouth,
                   pixels=None,
                   row_size=10,
                   edge_size=10,
                   strip_size=20,
                   create_fig=True,
                   show_fig=False,
                   return_imgdata=False):
    '''
    Find the normatlisation relative to the mean signal of CCD1 (red). The
    chip gaps are taken into account in order to find the appropriate
    response due to global flux variation, e.g. vignetting. The row_size,
    edge_size and strip_size are the number of pixels AFTER binning. The lower
    bounds are rounded down, upper bounds are rounded up.

    Parameters
    ----------
    flat: numpy.array
        The three reconstructed flat frames.
    binning: numeric
        The binning factor, which should be an integer, but data type is not
        important.
    northsout: str
        Indicate whether the Gemini north or south, because the chip gaps are
        different.
    pixels: numpy.array
        The pixel values adjusted for chip gaps for the corresponding pixel.
    row_size: numeric
        The number of rows on either side from the centre.
    edge_size: numeric
        The number of unusable edge columns.
    strip_size: numeric
        The number of columns to be used for fitting the straight line.
    create_fig: boolean
        Set to True to create a figure with matplotlib in the backend.
    show_fig: boolean
        Set to True to view the figure.
    return_imgdata: boolean
        Set to True to return an svg image as a BytesIO object.

    Returns
    -------
    flat: numpy.array
        The 3 individually normalised flats.
    normalisation: float
        The 2 ADU ratio: CCD2 to CCD1 and CCD3 to CCD2.
    imgdata: BytesIO
        BytesIO object holding the svg image.
    '''

    if northsouth == 'north':
        gap_binned_width = hamamatsu_n_gap / binx
    elif northsouth == 'south':
        gap_binned_width = hamamatsu_s_gap / binx
    else:
        raise ValueError('Please choose from "north" or "south".')

    if pixels == None:
        pixels = create_pixel_array(northsouth, binx)

    # Adjust for the binning
    binned_width = hamamatsu_width / binx
    row_middle = int(np.floor(hamamatsu_height / biny / 2))
    row_start = row_middle - int(row_size)
    row_end = row_middle + int(row_size)
    start_strip_start = int(edge_size)
    start_strip_end = int(np.ceil(start_strip_start + strip_size))
    end_strip_end = int(np.ceil(binned_width) - int(edge_size))
    end_strip_start = int(np.floor(end_strip_end - strip_size))
    gap = int(
        np.ceil(binned_width - end_strip_end + start_strip_start +
                gap_binned_width))

    # Normalise the individual flat frames
    flat1 = flat[0] / np.nanmean(flat[0])
    flat2 = flat[1] / np.nanmean(flat[1])
    flat3 = flat[2] / np.nanmean(flat[2])

    # Get the pixels near the chip gaps
    flat1_end_strip = np.mean(flat1[row_start:row_end],
                              axis=0)[end_strip_start:end_strip_end]
    flat2_start_strip = np.mean(flat2[row_start:row_end],
                                axis=0)[start_strip_start:start_strip_end]
    flat2_end_strip = np.mean(flat2[row_start:row_end],
                              axis=0)[end_strip_start:end_strip_end]
    flat3_start_strip = np.mean(flat3[row_start:row_end],
                                axis=0)[start_strip_start:start_strip_end]

    # Get a rough guess of the ADU ratio between the adjacent CCDs
    ratio_2to1 = np.mean(flat2_start_strip / flat1_end_strip)
    ratio_3to2 = np.mean(flat3_start_strip / flat2_end_strip)

    # Compute the normalisation by fitting a straight line
    normalisation2 = optimize.minimize(
        linear_fit, (-0.001, ratio_2to1, flat1_end_strip[0]),
        args=(flat1_end_strip, flat2_start_strip, gap),
        method='Nelder-Mead',
        options={
            'maxiter': 100000
        }).x[1]
    normalisation3 = optimize.minimize(
        linear_fit, (-0.001, ratio_3to2, flat2_end_strip[0]),
        args=(flat2_end_strip, flat3_start_strip, gap),
        method='Nelder-Mead',
        options={
            'maxiter': 100000
        }).x[1] * normalisation2

    # plot the normalised response function
    if create_fig:
        fig = plt.figure()
        plt.gcf().set_size_inches(10, 2.3)
        plt.clf()
        plt.plot(pixels[:int(len(pixels) / 3)],
                 flat1[row_middle],
                 color='red',
                 label='Red CCD')
        plt.plot(pixels[int(len(pixels) / 3):int(len(pixels) / 3 * 2)],
                 (flat2 / normalisation2)[row_middle],
                 color='green',
                 label='Green CCD')
        plt.plot(pixels[int(len(pixels) / 3 * 2):],
                 (flat3 / normalisation3)[row_middle],
                 color='blue',
                 label='Blue CCD')
        plt.xlim(min(pixels), max(pixels))
        plt.ylim(
            0.,
            np.nanpercentile(
                np.concatenate(
                    (flat1[row_middle], (flat2 / normalisation2)[row_middle],
                     (flat3 / normalisation3)[row_middle])), 99.9))
        plt.xlabel('Spectral Direction / Pixel')
        plt.ylabel('Relative Response')
        plt.grid()
        plt.tight_layout()

        if show_fig:
            plt.show()

        if return_imgdata:
            imgdata = BytesIO()
            fig.savefig(imgdata, format='svg')
            imgdata.seek(0)  # rewind the data
            plt.close()
            return (flat1, flat2, flat3), (normalisation2,
                                           normalisation3), imgdata
        else:
            plt.close()

    else:
        return (flat1, flat2, flat3), (normalisation2, normalisation3)


def flatten_image(image,
                  flat,
                  binx,
                  biny,
                  normalisation,
                  create_fig=True,
                  log=True,
                  show_gap=True,
                  show_fig=False,
                  return_imgdata=False):
    '''
    Field flattening the 3 images, correct for the response to the same level.

    Parameters
    ----------
    image: numpy.array
        The 3 reconstructed images for each of the CCDs.
    flat: numpy.array
        The 3 individually normalised flats.
    binning: numeric
        The binning factor, which should be an integer, but data type is not
        important.
    normalisation: float
        The 2 ADU ratio: CCD2 to CCD1 and CCD3 to CCD2.
    create_fig: boolean
        Set to True to create a figure with matplotlib in the backend.
    log: boolean
        Set to True to log the ADUs for better contrast.
    show_gap: boolean
        Set to True to display the chip gaps (arbitrary separation).
    show_fig: boolean
        Set to True to view the figure.
    return_imgdata: boolean
        Set to True to return an svg image as a BytesIO object.

    Returns
    -------
    out_image: numpy.array
        The flattened image.
    imgdata: BytesIO
        BytesIO object holding the svg image.
    '''

    flat1 = flat[0] / np.nanmean(flat[0])
    flat2 = flat[1] / np.nanmean(flat[1])
    flat3 = flat[2] / np.nanmean(flat[2])
    out_image = np.column_stack(
        (image[0] / flat1, image[1] / (flat2 / normalisation[0]), image[2] /
         (flat3 / normalisation[1])))[hamamatsu_middle_chip_min //
                                      biny:hamamatsu_middle_chip_max // binx]
    if create_fig:
        if show_gap:
            fig, axes = plt.subplots(1,
                                     3,
                                     figsize=(10, 4.5),
                                     sharex=True,
                                     sharey=True)
            plt.gcf().set_size_inches(10, 4.5)
            width = int(len(image[0]) / 3)
            if log:
                axes[0].imshow(np.log10(out_image[:, :width]),
                               aspect='auto',
                               origin='lower')
                axes[1].imshow(np.log10(out_image[:, width:width * 2]),
                               aspect='auto',
                               origin='lower')
                axes[2].imshow(np.log10(out_image[:, width * 2:width * 3]),
                               aspect='auto',
                               origin='lower')
            else:
                axes[0].imshow(out_image[:, :width],
                               aspect='auto',
                               origin='lower')
                axes[1].imshow(out_image[:, width:width * 2],
                               aspect='auto',
                               origin='lower')
                axes[2].imshow(out_image[:, width * 2:width * 3],
                               aspect='auto',
                               origin='lower')
            axes[1].set_title('Reconstructed Image')
            axes[1].set_xlabel('Spectral Direction / Pixel')
            axes[0].set_ylabel('Spatial Direction / Pixel')
        else:
            fig = plt.figure()
            plt.gcf().set_size_inches(10, 2.3)
            plt.imshow(np.log10(out_image), aspect='auto', origin='lower')
            plt.title('Reconstructed Image')
            plt.xlabel('Spectral Direction / Pixel')
            plt.ylabel('Spatial Direction / Pixel')
        plt.tight_layout()
        if show_fig:
            plt.show()

        if return_imgdata:
            imgdata = BytesIO()
            fig.savefig(imgdata, format='svg')
            imgdata.seek(0)  # rewind the data
            plt.close()
            return out_image, imgdata
        else:
            plt.close()
            return out_image
    else:
        return out_image
