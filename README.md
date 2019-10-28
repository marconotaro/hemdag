## Brief Description
**HEMDAG** package:

- implements several Hierarchical Ensemble Methods (HEMs) for Directed Acyclic Graphs (DAGs);
- reconciles flat predictions with the topology of the ontology;
- can enhance predictions of virtually any flat learning methods by taking into account the hierarchical relationships between ontology classes;
- provides biologically meaningful predictions that obey the _true-path-rule_, the biological and logical rule that governs the internal coherence of biomedical ontologies;
- is specifically designed for exploiting the hierarchical relationships of DAG-structured taxonomies, such as the Human Phenotype Ontology (HPO) or the Gene Ontology (GO), but can be safely applied to tree-structured taxonomies as well (as FunCat), since trees are DAGs;
- scales nicely both in terms of the complexity of the taxonomy and in the cardinality of the examples;
- provides several utility functions to process and analyze graphs;
- provides several performance metrics to evaluate HEMs algorithms.

## Documentation
Please get a look to the [documentation](http://hemdag-tutorials.readthedocs.io "HEMDAG’s documentation") to know how to download, install and make experiments with the **HEMDAG** package. 

## Cite HEMDAG
If you use HEMDAG, please cite our [BMC Bioinformatics article](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1854-y):

```
M. Notaro, M. Schubach, P. N. Robinson, and G Valentini. 
Prediction of Human Phenotype Ontology terms by means of hierarchical ensemble methods.
BMC Bioinformatics, 18(1):449, 2017
```
