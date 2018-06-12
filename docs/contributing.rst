.. _contributing:

============
Contributing
============

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

You can contribute in many ways:

----------------------
Types of Contributions
----------------------


Report Bugs
===========

Report bugs at https://github.com/visze/HEMDAG/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.


Fix Bugs
========

Look through the Github issues for bugs.
If you want to start working on a bug then please write short message on the issue tracker to prevent duplicate work.


Implement Features
==================

Look through the Github issues for features.
If you want to start working on an issue then please write short message on the issue tracker to prevent duplicate work.


Write Documentation
===================

HEMDAG could always use more documentation, whether as part of the official HEMDAG docs, in docstrings, or even on the web in blog posts, articles, and such.

HEMDAG uses `Sphinx <https://sphinx-doc.org>`_ for the user manual (that you are currently reading).
See `doc_guidelines` on how the documentation reStructuredText is used.
See `doc_setup` on creating a local setup for building the documentation.


Submit Feedback
===============

The best way to send feedback is to file an issue at https://github.com/visze/HEMDAG/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :)


.. _doc_guidelines:

------------------------
Documentation Guidelines
------------------------

For the documentation, please adhere to the following guidelines:

- Put each sentence on its own line, this makes tracking changes through Git SCM easier.
- Provide hyperlink targets, at least for the first two section levels.
- Use the section structure from below.

::

    .. heading_1:

    =========
    Heading 1
    =========


    .. heading_2:

    ---------
    Heading 2
    ---------


    .. heading_3:

    Heading 3
    =========


    .. heading_4:

    Heading 4
    ---------


    .. heading_5:

    Heading 5
    ~~~~~~~~~


    .. heading_6:

    Heading 6
    :::::::::


.. _doc_setup:

-------------------
Documentation Setup
-------------------

For building the documentation, you have to install the Python program Sphinx.
We use conda for that, see :ref:`conda`

Use the following steps for installing Sphinx and the dependencies for building the HEMDAG documentation:

.. code-block:: console

    $ cd HEMDAG/docs
    $ conda create --name sphinx --file environment.yml
    $ source activate sphinx

Use the following for building the documentation.
If you are not in the sphinx environment (e.g. you uses ``source deactivate sphinx``) please activate the virtual environment using ``source activate sphinx``
Afterwards, you can always use ``make html`` for building.

.. code-block:: console

    (sphinx) $ cd HEMDAG/docs
    (sphinx) $ make html  # rebuild for changed files only
    (sphinx) $ make clean && make html  # force rebuild


------------
Get Started!
------------

Ready to contribute?

1. Fork the `HEMDAG` repo on GitHub.
2. Clone your fork locally:

  .. code-block:: console

      $ git clone git@github.com:your_name_here/HEMDAG.git

3. Create a branch for local development:

  .. code-block:: console

      $ git checkout -b name-of-your-bugfix-or-feature

  Now you can make your changes locally.

5. When you're done making your changes, make sure that the build runs through.
   For R:

   .. code-block:: console

    $ TODO

   For documentation:

   .. code-block:: console

      $ cd docs && make clean && make html

6. Commit your changes and push your branch to GitHub:

  .. code-block:: console

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.


-----------------------
Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated.
3. Describe your changes in the ``CHANGELOG`` file.
