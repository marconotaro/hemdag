.. _quickstart:

============
Quickstart
============

This short *HowTo* guides you from downloading HEMDAG library, load it into your R environment and make first computations.

Installation
================

Please go to the :ref:`installation` section and install HEMDAG by using one of the ways shown.

Load HEMDAG Library
=======================

Start R in your console using

.. code-block:: console

    $ R

then load the library by using

.. code-block:: R

    > library(HEMDAG)

First Classification -- for the Impatient
=============================================

HEDMDAG encompasses in total 23 hierarchical ensemble methods. Below we show the *simple* call to all the hierarchical ensemble algorithms included in HEMDAG, bu using the pre-built datasets available in the HEMDAG for making predictions. For more details about datasets and methods have a look to section :ref:`tutorial`.

A. Loading the pre-built dataset of HEMDAG

.. code-block:: R

    # load the ontology DAG g
    > data(graph);

    # load the scores matrix S
    > data(scores);

    # load the annotation matrix L
    > data(labels);

    # compute the root node
    > root <- root.node(g);


B. HTD-DAG: Hierarchical Top-Down for DAG

.. code-block:: R

    > S.htd  <- htd(S, g, root);


C. GPAV: Generalized Pool-Adjacent-Violators

.. code-block:: R

    > S.gpav <- gpav.over.examples(S, g, W=NULL);


D. TPR-DAG (True Path Rule for DAG) and all its 18 ensemble variants

.. code-block:: R

    > S.tprTF         <- tpr.dag(S, g, root, positive="children", bottomup="threshold.free", topdown="htd");
    > S.tprT          <- tpr.dag(S, g, root, positive="children", bottomup="threshold", topdown="htd", t=0.5);
    > S.tprW          <- tpr.dag(S, g, root, positive="children", bottomup="weighted.threshold.free", topdown="htd", w=0.5);
    > S.tprWT         <- tpr.dag(S, g, root, positive="children", bottomup="weighted.threshold", topdown="htd", t=0.5, w=0.5);

    > S.descensTF     <- tpr.dag(S, g, root, positive="descendants", bottomup="threshold.free", topdown="htd");
    > S.descensT      <- tpr.dag(S, g, root, positive="descendants", bottomup="threshold", topdown="htd", t=0.5);
    > S.descensW      <- tpr.dag(S, g, root, positive="descendants", bottomup="weighted.threshold.free", topdown="htd", w=0.5);
    > S.descensWT     <- tpr.dag(S, g, root, positive="descendants", bottomup="weighted.threshold", topdown="htd", t=0.5, w=05);
    > S.descensTAU    <- tpr.dag(S, g, root, positive="descendants", bottomup="tau", topdown="htd", t=0.5);

    > S.isotprTF      <- tpr.dag(S, g, root, positive="children", bottomup="threshold.free", topdown="gpav");
    > S.isotprT       <- tpr.dag(S, g, root, positive="children", bottomup="threshold", topdown="gpav", t=0.5);
    > S.isotprW       <- tpr.dag(S, g, root, positive="children", bottomup="weighted.threshold.free", topdown="gpav", w=0.5);
    > S.isotprWT      <- tpr.dag(S, g, root, positive="children", bottomup="weighted.threshold", topdown="gpav", t=0.5, w=0.5);

    > S.isodescensTF  <- tpr.dag(S, g, root, positive="descendants", bottomup="threshold.free", topdown="gpav");
    > S.isodescensT   <- tpr.dag(S, g, root, positive="descendants", bottomup="threshold", topdown="gpav", t=0.5);
    > S.isodescensW   <- tpr.dag(S, g, root, positive="descendants", bottomup="weighted.threshold.free", topdown="gpav", w=0.5);
    > S.isodescensWT  <- tpr.dag(S, g, root, positive="descendants", bottomup="weighted.threshold", topdown="gpav", t=0.5, w=0.5);
    > S.isodescensTAU <- tpr.dag(S, g, root, positive="descendants", bottomup="tau", topdown="gpav", t=0.5);


E. Obozisnki heuristic methods

.. code-block:: R

    > S.max <- obozinski.max(S,g,root);
    > S.and <- obozinski.and(S,g,root);
    > S.or  <- obozinski.or(S,g,root);
