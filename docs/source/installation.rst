.. BHTOM-SPEC documentation master file, created by
   sphinx-quickstart on Thu Jan  9 11:38:33 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation Guide
==================

Dependencies
^^^^^^^^^^^^

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

    You can install by cloning and then execute the installation

    .. code-block:: shell

        git clone https://github.com/cylammarco/bhtomspec
        pip3 install -e bhtomspec

    Or install directly from the github respository

    .. code-block:: shell

        pip3 install https://github.com/cylammarco/bhtomspec/archive/master.zip

    You may need to use ``sudo`` or the ``--user`` flag.

    In some cases, it may run into the error ``Could not install packages due to an EcvirnomentalError: [Errno 2]``. Try running the pip command again, the problem usually goes away. It may also help by including the upgrade and fore-reinstall flag: ``pip3 install --upgrade --force-reinstall -e bhtomspec``.

Other repositories for BHTOM Spec
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is not currently registered on any package index service. It may happen if this project carries on to a full integration with BHTOM.
