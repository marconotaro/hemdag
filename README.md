# Welcome to HEMDAG R package!
|          | Bioconda Badges | CRAN  Badges |
| -------- | --------------- | ------------ |
| License | [![GPL3-License](https://img.shields.io/cran/l/HEMDAG?label=License&style=flat-square&color=orange)](https://www.gnu.org/licenses/gpl-3.0.en.html) | [![GPL3-License](https://img.shields.io/cran/l/HEMDAG?label=License&color=orange)](https://www.gnu.org/licenses/gpl-3.0.en.html)
| Last Version | [![Anaconda-Version](https://anaconda.org/bioconda/r-hemdag/badges/version.svg)](https://anaconda.org/bioconda/r-hemdag) | [![CRAN-version](http://www.r-pkg.org/badges/version/HEMDAG?color=blue)](https://cran.r-project.org/package=HEMDAG)
| Last Updated | [![Anaconda-Last-Updated](https://anaconda.org/bioconda/r-hemdag/badges/latest_release_date.svg)](https://anaconda.org/bioconda/r-hemdag)| [![CRAN-Last-Updated](https://www.r-pkg.org/badges/last-release/HEMDAG?color=blue)](https://cran.r-project.org/package=HEMDAG)
| Read the Docs | [![Read the Docs Status](https://img.shields.io/readthedocs/hemdag?logo=read%20the%20docs&logoColor=white&style=flat-square)](https://hemdag.readthedocs.io) | [![Read the Docs Status](https://img.shields.io/readthedocs/hemdag?logo=read%20the%20docs&logoColor=white)](https://hemdag.readthedocs.io)
| Total Downloads | [![Anaconda-Downloads](https://anaconda.org/bioconda/r-hemdag/badges/downloads.svg)](https://anaconda.org/bioconda/r-hemdag) | [![CRAN-Downloads](http://cranlogs.r-pkg.org/badges/grand-total/HEMDAG?color=green)](https://cranlogs.r-pkg.org/downloads/total/2017-08-11:last-day/HEMDAG)
| Code Coverage | [![CodeCov](https://img.shields.io/codecov/c/gh/marconotaro/hemdag?logo=codecov&style=flat-square)](https://codecov.io/gh/marconotaro/hemdag) | [![CodeCov](https://img.shields.io/codecov/c/gh/marconotaro/hemdag?logo=codecov)](https://codecov.io/gh/marconotaro/hemdag)
<!-- trick: to know daily HEMDAG downloads (and from how many days HEMDAG lives) replace total with daily in the cranlogs link -->

---

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
Please get a look to the [documentation](https://hemdag.readthedocs.io "HEMDAG’s documentation") to know how to download, install and make experiments with the **HEMDAG** package.
## Cite HEMDAG
If you use HEMDAG, please cite our [Bioinformatics article](https://doi.org/10.1093/bioinformatics/btab485) or [BMC Bioinformatics article](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1854-y):

```
Marco Notaro, Marco Frasca, Alessandro Petrini, Jessica Gliozzo, Elena Casiraghi, Peter N Robinson, Giorgio Valentini
HEMDAG: a family of modular and scalable hierarchical ensemble methods to improve Gene Ontology term prediction,
Bioinformatics, Volume 37, Issue 23, 1 December 2021, Pages 4526–4533

M. Notaro, M. Schubach, P. N. Robinson, and G Valentini.
Prediction of Human Phenotype Ontology terms by means of hierarchical ensemble methods.
BMC Bioinformatics, 18(1):449, 2017
```
