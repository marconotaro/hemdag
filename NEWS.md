#### HEMDAG 2.2.5

###### NEW FEATURES
- precision-recall performance computed through ``precrec`` package:
	- added ``precision.at.all.recall.levels.single.class``;
	- ``PXR.at.multiple.recall.levels.over.classes`` substituted with ``precision.at.given.recall.levels.over.classes``;
- improved IO functions: the extension of the input or output file can be or plain text (``.txt``) or compressed (``.gz``); 

###### CHANGES
- fixed minor bugs;
- manual improved;

#### HEMDAG 2.2.4

###### CHANGES
- fixed ``CRAN`` Package Check Results: removed unneeded header and define from ``GPAV C++`` source code

#### HEMDAG 2.2.3

###### NEW FEATURES
- Added ``GPAV`` algorithm (Burdakov et al., *Journal of Computational Mathematics*, 2006 -- [link](https://www.jstor.org/stable/43693336));
- Embedded ``GPAV`` algorithm in the top-down step of the functions ``TPR.DAG``, ``Do.TPR.DAG`` and ``Do.TPR.DAG.holdout``;
- Some functions have been defunct. To know the defunct functions just typing in the R environment: ``help("HEMDAG-defunct")``;

###### CHANGES
- manual improved;

###### AUTHOR
- Added **Alessandro Petrini** as author for his contribution in writing the ``C++`` code of ``GPAV`` algorithm;

#### HEMDAG 2.1.3

###### CHANGES
- various fixes from 2.1.2

#### HEMDAG 2.1.2

###### NEW FEATURES
- Improved performance metrics:
	- added ``compute.Fmeasure.multilabel``;
	- added ``PXR.at.multiple.recall.levels.over.classes``;
	- all the performance metrics (``AUPRC``, ``AUROC``, ``FMM``, ``PXR``) can be computed either **one-shot** or averaged **across folds**;

- Improved the high-level hierarchical ensemble functions:
	- embedded the new performance metric functions;
	- added the parameter ``metric``: maximization by ``FMAX`` or ``PRC`` (see manual for further details);
	- added some checkers (warning/stop messages) to make the library more user-friendly;

###### CHANGES
- manual improved;

#### HEMDAG 2.0.1

###### CHANGES
- corrected bug in ``do.stratified.cv.data.single.class``;

#### HEMDAG 2.0.0

###### NEW FEATURES
- Added ``TPR-DAG``: function gathering several hierarchical ensemble variants;
- Added ``Do.TPR.DAG``: high-level function to run ``TPR-DAG`` **cross-validated** experiments;
- Added ``Do.TPR.DAG.holdout``: high-level functions to run ``TPR-DAG`` **holdout** experiments;

- The following ``TPR-DAG`` and ``DESCENS`` high-level functions were removed:
	- Do.tpr.threshold.free;
	- Do.tpr.threshold.cv;
	- Do.tpr.weighted.threshold.free.cv;
	- Do.tpr.weighted.threshold.cv;
	- Do.descens.threshold.free;
	- Do.descens.threshold.cv;
	- Do.descens.weighted.threshold.free.cv;
	- Do.descens.tau.cv;
	- Do.descens.weighted.threshold.cv;
	- Do.tpr.threshold.free.holdout;
	- Do.tpr.threshold.holdout;
	- Do.tpr.weighted.threshold.free.holdout;
	- Do.tpr.weighted.threshold.holdout;
	- Do.descens.threshold.free.holdout;
	- Do.descens.threshold.holdout;
	- Do.descens.weighted.threshold.free.holdout;
	- Do.descens.tau.holdout;
	- Do.descens.weighted.threshold.holdout;

> NOTE: all the removed functions can be run opportunely setting the input parameters of the new high-level function ``Do.TPR.DAG`` (for **cross-validated** experiments) and ``Do.TPR.DAG.holdout`` (for **hold-out** experiments);

###### CHANGES
- manual improved;
	
#### HEMDAG 1.1.1

###### NEW FEATURES
- Added ``DESCENS`` algorithm;
- Added Heuristic Methods ``MAX``, ``AND``, ``OR`` (Obozinski et al., Genome Biology, 2008 -- [link](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-s1-s6));
- Added ``tupla.matrix`` function;

###### CHANGES
- manual improved;
- Added link to the GitHub repository **HPOparser** ([link](https://github.com/marconotaro/HPOparser));
- Added ``CITATION`` file;

#### HEMDAG 1.0.0

###### PACKAGE GENESIS
