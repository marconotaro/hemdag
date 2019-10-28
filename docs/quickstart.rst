.. _quickstart:

============
Quickstart
============

This short How-To guides you from downloading the HEMDAG, load it into your R environment and make a first computation.

Installation
=====================

Please goto the :ref:`installation` section and use the :ref:`conda` option to install HEMDAG.

Load HEMDAG library
==========================

Start R in your console using

.. code-block:: console

    $ R

then load the library by using

.. code-block:: R

    library("HEMDAG")

Your first classification
==========================

We will use the DESCENS algorithm to do some predictions on a DAG (Human Phenotype Ontology).

In contrast to the *vanilla* TPR-DAG version, DESCENS takes into account the contribution of all the descendants of each node instead of only that of its children. So DESCENS predictions are more influenced by the information embedded in the most specific terms of the taxonomy (e.g. leaf nodes), thus putting more emphasis on the terms that most characterize the gene under study.

.. code-block:: R

    # load a ontology DAG stored in g
    data(graph);
    # load scores for genes to HPO stored in S
    data(scores);
    # load labels in L (genes annotated to HPO terms)
    data(labels);
    # Set the root in your DAG
    root <- root.node(g);
    # Run DESCENS Threshold free
    S.descensTF  <- TPR.DAG(S, g, root, positive="descendants", bottomup="threshold.free", topdown="HTD");
    # Run DESCENS with threshold
    S.descensT   <- TPR.DAG(S, g, root, positive="descendants", bottomup="threshold", topdown="HTD", t=0.5);
    # Run weighted DESCENS with threshold free
    S.descensW   <- TPR.DAG(S, g, root, positive="descendants", bottomup="weighted.threshold.free", topdown="HTD", w=0.5);
    # Run weighted DESCENS with threshold
    S.descensWT  <- TPR.DAG(S, g, root, positive="descendants", bottomup="weighted.threshold", topdown="HTD", t=0.5, w=05);
    # Run DESCENS TAU
    S.descensTAU <- TPR.DAG(S, g, root, positive="descendants", bottomup="tau", topdown="HTD", t=0.5);
