.. HEMDAG documentation master file, created by
   sphinx-quickstart on Fri May 11 12:02:41 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to HEMDAG's documentation!
====================================

HEMDAG package:

  * implements several Hierarchical Ensemble Methods (HEMs) for Directed Acyclic Graphs (DAGs);
  * reconciles flat predictions with the topology of the ontology;
  * can enhance the predictions of virtually any flat learning methods by taking into account the hierarchical relationships between ontology classes;
  * guarantees biologically meaningful predictions that always obey the *true-path-rule*, the biological and logical rule that governs the internal coherence of biomedical ontologies;
  * is specifically designed for exploiting the hierarchical relationships of DAG-structured taxonomies, such as the Human Phenotype Ontology (HPO) or the Gene Ontology (GO), but can be safely applied to tree-structured taxonomies as well (e.g. FunCat), since trees are DAGs;
  * scales nicely both in terms of the complexity of the taxonomy and in the cardinality of the examples;
  * provides several utility functions to process and analyze graphs;
  * provides several performance metrics to evaluate HEMs algorithms.

.. toctree::
   :caption: Installation & Getting Started
   :name: getting-started
   :maxdepth: 1
   :hidden:

   quickstart
   install

.. toctree::
    :caption: Usage & Tutorial
    :name: HEMDAG-usage
    :maxdepth: 1
    :hidden:

    usage
    tutorial

.. toctree::
    :caption: Tips & Tricks
    :name: tips-tricks
    :maxdepth: 1
    :hidden:

    faq


.. toctree::
    :caption: Project Info
    :name: project-info
    :maxdepth: 1
    :hidden:

    citing
    contributing
    authors
    history
    license
