.. _installation:

============
Installation
============

HEMDAG is available on CRAN as well as through Bioconda and also from source code. You can use one of the following ways for installing HEMDAG.


.. _conda:

Installation via Conda
========================

.. note::

    This is the recommended way of installing for normal users.

This is the recommended way to install HEMDAG because it will enables you to switch software versions easily. And in addition R with all needed dependencies will be installed.

First, you have to install the Miniconda Python3 distribution. See `here <https://conda.io/docs/install/quick.html>`_ for installation instructions. Make sure to ...

 - Install the *Python 3* version of Miniconda.
 - Answer yes to the question whether conda shall be put into your PATH.

Then, you can install HEMDAG with

.. code-block:: console

    $ conda install -c bioconda -c conda-forge r-hemdag

from the `Bioconda <https://bioconda.github.io>`_ channel.

Global Installation
========================

You can directly install the library via R by issuing

.. code-block:: console

    $ R -e "install.packages('HEMDAG',repos = 'http://cran.us.r-project.org')"

in your terminal. But be sure to install R properly before that command.

Alternatively, you can install the ``HEMDAG`` library by typing in the ``R`` environment:

 .. code-block:: R

	install.packages("HEMDAG");

Another possibility to install the development version of ``HEMDAG`` is by using the ``devtools`` package:

 .. code-block:: R

	library(devtools);
	install_github("gecko515/HEMDAG");

Before running the above commands be sure to have correctly installed the ``devtools`` package (`link <https://cran.r-project.org/web/packages/devtools/README.html>`_)


.. _install_from_source:

Installing from Source
=======================


This section describes how to build HEMDAG from scratch.

Prequisites
-----------

For building HEMDAG, you will need the following dependencies

 - R (≥ 2.10)
 - R-libraries:
    - PerfMeas
    - rbgl (bioconductor)
    - graph (bioconductor)
    - precrec
    - preprocessCore  (bioconductor)

Package from CRAN
-----------------------------------
On a linux environment, download the package source from the `CRAN repo <https://cran.rstudio.com/web/packages/HEMDAG/>`_ and save it in the folder ``pippo``. Then type:

 .. code-block:: console

	R CMD INSTALL pippo/HEMDAG_2.1.3.tar.gz



Direct Git Checkout
--------------------

.. note::

    You only need to install from source if you want to develop HEMDAG yourself.

In this tutorial, we will download the HEMDAG sources and build them in ``~/HEMDAG``:

.. code-block:: console

  ~ $ cd ~
  ~ $ git clone https://github.com/gecko515/HEMDAG.git HEMDAG

Building
--------

You can build HEMDAG using::

.. code-block:: console

  R CMD build HEMDAG

This will generate the file ``HEMDAG_2.1.3.tar.gz`` and just install the package via:

.. code-block:: console

  R CMD INSTALL HEMDAG_2.1.3.tar.gz
