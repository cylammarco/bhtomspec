.. _reduction_config:

Reduction Configuration
=======================

This only applies to GMOS reduction where a full reduction is performd. SPRAT reductions begin with flattened spectral image data.

GMOS field-flattening configuration (YAML)
------------------------------------------

.. code-block:: yaml
   :linenos:

   folder_path: "."
   output_folder_path: "output"
   light_folder: ~
   flat_folder: ~
   arc_folder: ~
   bias_folder: ~
   arc_flat_folder: ~
   arc_bias_folder: ~
   light_extension: ["fits", "fit", "bz2", "gz"]
   flat_extension: ["fits", "fit", "bz2", "gz"]
   arc_extension: ["fits", "fit", "bz2", "gz"]
   bias_extension: ["fits", "fit", "bz2", "gz"]
   arc_flat_extension: ["fits", "fit", "bz2", "gz"]
   arc_bias_extension: ["fits", "fit", "bz2", "gz"]
   bias_master_filename: bias_master
   arc_bias_master_filename: arc_bias_master
   flat_row_size: 10
   flat_edge_size: 10
   flat_strip_size: 20
   force_recreate_bias: false
   save_bias: true
   save_bias_format: ["fits", "npy"]
   overwrite_bias_image: true
   save_flattened_image: true
   save_flattened_image_format: ["fits", "npy"]
   overwrite_flattened_image: true
   create_fig: true
   show_fig: false
   diagnostic_pdf: true

The `light_folder`, `flat_folder`, `arc_folder`, `bias_folder`, `arc_flat_folder` and `arc_bias_folder` are relative to the `folder_path`. If they are empty, the folders are assumed to be called `light`, `flat`, `arc`, `bias`, `arc_flat`, `arc_bias` folders located at `folder_path`.

All the files with the extension listed in the respective extension arguments (`light_extension`...) will be processed.

`flat_row_size` is the number of rows on either size from the central row of the spectral image that are to be summed to compute the sensitivity function. `flat_edge_size` is the number of edge pixels that are removed in the computation. `flat_strip_size` is the number of pixels from the `flat_edge_size`-th pixel to be fitted a straight line over.
