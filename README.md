# Black Hole TOM Spectroscopic Reduction

[![Documentation Status](https://readthedocs.org/projects/bhtom-spec/badge/?version=latest)](https://bhtom-spec.readthedocs.io/en/latest/?badge=latest)

This provides the spectroscopic data reduction for BHTOM. At the time of
writing, it is planned for implementing the longslit spectral reduction of
Gemini/GMOS and LT/SPRAT.

## Dependencies

_Image fieldflattening_
* astropy
* numpy
* pyyaml
* reportlab (Optional for pdf diagnostic plots)
* svglib (Optional for pdf diagnostic plots)
* scipy

_Spectral extraction and calibration_ (Nor required if to perform fieldflattening only)
* aspired >= 0.0.1 (Spectral Extraction)
* rascal >= 0.1 (For wavelength calibration)
* spectres >= 2.1.1

At the time of writing, the aspired and rascal should be installed from the respective development branches
* https://github.com/cylammarco/aspired.git@dev
* https://github.com/jveitchmichaelis/rascal.git@dev

This is also required if you wish to build the document locally

+ https://github.com/SuperKogito/sphinxcontrib-pdfembed

## GMOS Longslit Reduction

This currently allows semi-automated reduction for the longslit mode of both GMOS North and South. See more about GMOS [here](https://www.gemini.edu/instrumentation/current-instruments/gmos).

- [x] Reconstruct the GMOS Hamamatsu fits into a 2D image
- [x] Create master bias
- [x] Normalise response function of the flats from 3 CCDs
- [x] Flatten the image
- [x] Handles binning automatically
- [x] Handles _mismatched_ binning between flattened 2D spectral data image and 2D arc image *only*
- [x] Handle 2 types of ROI: _Full Frame_ and _Central Spectrum_

N.B. The central spectrum ROIs are reconstructed into full frames.

## SPRAT Longslit Reduction

See more about SPRAT [here](http://telescope.livjm.ac.uk/TelInst/Inst/SPRAT/)

- [x] Re-reduction based only on the final multiple extension data product
- [x] Re-reduction based on the final multiple extension data product and an arc image
- [x] Generic spectral reduction
