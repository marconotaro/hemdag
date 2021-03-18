.. _installation:

============
Installation
============
HEMDAG is available on CRAN, through Bioconda and from source code. You can use one of the following ways for installing HEMDAG.

.. _conda:

Installation via Conda
========================
This is the recommended way to install HEMDAG for normal users because it will enable you to switch software versions easily. In addition R with all needed dependencies will be installed.

First, you have to install the Miniconda Python3 distribution. See `here <https://docs.conda.io/en/latest/miniconda.html>`_ for installation instructions. Make sure to:

 - install the Python3 version of Miniconda.
 - answer yes to the question whether conda shall be put into your PATH.

Then, you can install HEMDAG with

.. code-block:: console

    $ conda install -c bioconda -c conda-forge r-hemdag

from the `Bioconda <https://bioconda.github.io>`_ channel.

Global Installation
========================
You can directly install the library via R by typing in your terminal:

.. code-block:: console

    $ R -e 'install.packages("HEMDAG", repos="http://cran.us.r-project.org")'

Alternatively, you can install the HEMDAG library by typing in the R environment:

.. code-block:: R

    install.packages("HEMDAG");

Another possibility to install the development version of HEMDAG is by using the *devtools* package (`link <https://CRAN.R-project.org/package=devtools>`_):

.. code-block:: R

    library(devtools);
    install_github("marconotaro/hemdag");

Dependencies
=======================
To install or build HEMDAG the following dependencies are required:

 - R (≥ 2.10)

 - R dependencies

    - graph (bioconductor)
    - rbgl (bioconductor)
    - precrec
    - preprocessCore  (bioconductor)
    - plyr
    - foreach
    - doParallel

.. note::

    CRAN does not automatically install Bioconductor packages. To install them:

    .. code-block:: R

        if(!requireNamespace("BiocManager", quietly=TRUE))
            install.packages("BiocManager")

        BiocManager::install("graph")
        BiocManager::install("rbgl")
        BiocManager::install("preprocessCore")

Installing from Source
=======================
Here we describe how to build HEMDAG from scratch.

Package from CRAN
-----------------------------------
On a Linux environment, download the package source from the `CRAN repo <https://CRAN.R-project.org/package=HEMDAG>`_ and save it (for instance) in the folder ``pippo``. Then type:

 .. code-block:: console

    R CMD INSTALL pippo/HEMDAG_<pkg-version-number>.tar.gz

.. note::

    Replace ``<pkg-version-number>`` with the version number of the downloaded HEMDAG package.

Direct Git Checkout
--------------------
.. note::

    You only need to install from source if you want to develop HEMDAG yourself.

Below, we will download the HEMDAG sources and build them in ``~/hemdag``:

.. code-block:: console

  ~ $ cd ~
  ~ $ git clone https://github.com/marconotaro/hemdag.git

Building
--------
You can build HEMDAG by using:

.. code-block:: console

  R CMD build hemdag

This will generate the file ``HEMDAG_<package-version-number>.tar.gz`` and just install the package via:

.. code-block:: console

  R CMD INSTALL HEMDAG_<package-version-number>.tar.gz
