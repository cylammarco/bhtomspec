.. BHTOM-SPEC documentation master file, created by
   sphinx-quickstart on Thu Jan  9 11:38:33 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation Guide
==================

Dependencies
^^^^^^^^^^^^

    Python >= 3.6

    These packages can be installed with pip

      + numpy >= 1.16
      + scipy >= 1.5
      + spectres >= 2.1.1
      + reportlab
      + svglib
      + pyyaml

    These have to be installed from the respective development branches at

      + https://github.com/cylammarco/aspired#dev
      + https://github.com/jveitchmichaelis/rascal#dev

    This is required if you wish to build the document locally

      + https://github.com/SuperKogito/sphinxcontrib-pdfembed

BHTOM Spec
^^^^^^^^^^

    Clone the repositry

    .. code-block:: shell

        git clone https://github.com/cylammarco/bhtomspec
        git clone https://github.com/cylammarco/bhtomspec-example

    To test the bhtomspec, all the examplex should run, e.g.

    .. code-block:: shell

        python3 [some_path_a]/bhtomspec/SPRAT/sprat_extraction.py [some_path_b]/bhtomspec-example/SPRAT/example/ExtractionCase1/20180810_lhs6328_case1.yaml
        python3 [some_path_a]/bhtomspec/GMOS/gmos_reduction.py [some_path_b]/bhtomspec-example/GMOS/example/ReductionCase1/flattening_config.yaml
        python3 [some_path_a]/bhtomspec/GMOS/gmos_extraction.py [some_path_b]/bhtomspec-example/GMOS/example/ExtractionCase1/extraction_config_1.yaml

    There are 3 shell scripts that run all the examples, however, the path has to be configured before they work. For example, in ``[some_path_b]/bhtomspec-example/SPRAT/run_all_examples.sh``

    .. code-block:: shell

        python3 sprat_extraction.py example/ExtractionCase1/20180810_lhs6328_case1.yaml
        python3 sprat_extraction.py example/ExtractionCase2/20180810_lhs6328_case2.yaml
        python3 sprat_extraction.py example/ExtractionCase3/20180810_lhs6328_case3.yaml
        python3 sprat_extraction.py example/ExtractionCase4/20180810_lhs6328_case4.yaml
        python3 sprat_extraction.py example/ExtractionCase5/20180810_lhs6328_case5.yaml

    has to be modified to

    .. code-block:: shell

        python3 [some_path_a]/bhtomspec/SPRAT/ssprat_extraction.py example/ExtractionCase1/20180810_lhs6328_case1.yaml
        python3 [some_path_a]/bhtomspec/SPRAT/ssprat_extraction.py example/ExtractionCase2/20180810_lhs6328_case2.yaml
        python3 [some_path_a]/bhtomspec/SPRAT/ssprat_extraction.py example/ExtractionCase3/20180810_lhs6328_case3.yaml
        python3 [some_path_a]/bhtomspec/SPRAT/ssprat_extraction.py example/ExtractionCase4/20180810_lhs6328_case4.yaml
        python3 [some_path_a]/bhtomspec/SPRAT/ssprat_extraction.py example/ExtractionCase5/20180810_lhs6328_case5.yaml

    before it can run.


Other repositories for BHTOM Spec
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is currently not a pacakge, just scripts, so it is not registered on any package index service. It may happen if this project continues onto a full integration with BHTOM.
