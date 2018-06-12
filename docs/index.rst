.. HEMDAG documentation master file, created by
   sphinx-quickstart on Fri May 11 12:02:41 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to HEMDAG's documentation!
==================================


*HEMDAG* library:

- implements several Hierarchical Ensemble Methods (HEM) for Directed Acyclic Graphs (DAGs);
- can be used to *reconcile* flat predictions by considering the topology of the ontology;
- provides biologically consistent predictions according to the hierarchical nature of the ontology classes;
- is specifically designed for exploiting the hierarchical relationships of DAG-structured taxonomies, such as the Human Phenotype Ontology (``HPO``) or the Gene Ontology (``GO``). However ``HEMDAG`` can be also safely applied to tree-structured taxonomies (as ``FunCat``), since trees are DAGs;
- scales nicely both in terms of the complexity of the taxonomy and in the cardinality of the examples.

.. toctree::
   :caption: Installation & Getting Started
   :name: getting-started
   :maxdepth: 1
   :hidden:

   quickstart
   install

.. toctree::
    :caption: Usage & Tutorials
    :name: HEMDAG-usage
    :maxdepth: 1
    :hidden:

    usage
    tutorials

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




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
