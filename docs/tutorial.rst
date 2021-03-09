.. role:: R(code)
   :language: R

.. _tutorial:

================================
Tutorial
================================
In this tutorial we show a step-by-step application of HEMDAG to the hierarchical prediction of associations between human gene and abnormal phenotype. To this end we will use the small pre-built dataset available in the HEMDAG library. However, all the hierarchical ensemble methods encompassed in HEMDAG library can be run by using:

    * any ontology listed in OBO foundry (`link <http://www.obofoundry.org>`__);
    * any flat score matrix, achieved by using any flat classifier ranging from linear, to probabilistic methods, to neural networks, to gradient boosting and many others;
    * any annotation matrix.

Of course, the number of terms among the graph, the flat scores matrix and the annotation matrix must match.

.. note::

    For the experiments shown below, we used the latest version of HEMDAG package, the R version 3.6.3 and on a machine having Ubuntu 18.04 as operative system.

Load the HEMDAG Library
==============================
To load the HEMDAG library, open the R environment and type:

.. code-block:: R

    > library(HEMDAG);

Load the Flat Scores Matrix
================================
In their more general form, the hierarchical ensemble methods adopt a two-step learning strategy:

    * the first step consists in the flat learning of the ontology terms;
    * the second step *reconciles* the flat predictions by considering the topology of the underlying ontology.

Consequently, the first *ingredient* that we need in a hierarchical ensemble classification is the flat scores matrix. For the sake of simplicity, in the examples shown below we use the pre-built dataset available in the HEMDAG library. To load the flat scores matrix, type in the the R environment:

.. code-block:: R

    > data(scores);

With the above command we loaded the flat scores matrix ``S``, that is a named 100 X 23 matrix. Rows correspond to genes (Entrez GeneID) and columns to HPO terms/classes. The scores represent the likelihood that a given gene belongs to a given class: the higher the value, the higher the likelihood that a gene belongs to a given class. This flat scores matrix was obtained by running the RANKS package (`link <https://cran.rstudio.com/web/packages/RANKS/>`__).

Normalization
----------------
Since RANKS **returns a score and not a probability**, we must normalize the scores of the matrix ``S`` to make the flat scores comparable with the hierarchical ones. In the case the flat classifier returns directly a probability there is no needed to normalize the flat scores matrix, since the flat scores can be directly compared with the hierarchical ones.

HEMDAG allows to normalize the flat scores according to two different procedures:

1. **maxnorm**: Normalization in the sense of the maximum: the score of each class is normalized by dividing the score values for the maximum score of that class:

.. code-block:: R

    > S.maxnorm <- scores.normalization(norm.type="maxnorm", S);

2. **qnorm**: Quantile normalization: quantile normalization of the *preprocessCore* package is used:

.. code-block:: R

    > S.qnorm <- scores.normalization(norm.type="qnorm", S);

Be sure to install the *preprocessCore* package before running the above command. Yo can install it by conda (``conda install -c bioconda bioconductor-preprocesscore``) or by Bioconductor (`link <https://bioconductor.org/packages/release/bioc/html/preprocessCore.html>`_)

.. note::

    For the sake of simplicity, in all the examples shown in section :ref:`hem`, the input flat scores matrix was normalized according to the ``maxnorm`` normalization:

    .. code-block:: R

        S.norm <- scores.normalization(norm.type="maxnorm", S);

Load the Graph
=================
In order to know how the hierarchical structure of the HPO terms, we need to load the graph:

.. code-block:: R

    > data(graph);

With the above command we loaded the graph ``g``, an object of class ``graphNEL``. The graph ``g`` has 23 nodes and 30 edges and represents the *ancestors view* of the HPO term ``Camptodactyly of finger`` (`HP:0100490 <https://hpo.jax.org/app/browse/term/HP:0100490>`_). Nodes of the graph ``g`` correspond to terms of the flat scores matrix ``S``.

Plot the Graph (optional)
-----------------------------
.. note::
    To plot the graph you need to install before the `Rgraphviz` package. Yo can install this library for example by conda (``conda install -c bioconda bioconductor-rgraphviz``) or by Bioconductor (`link <https://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html>`__).

If you want to visualize the *ancestors view* of the term ``HP:0100490``, just type:

.. code-block:: R

    > library(Rgraphviz);
    > plot(g);

.. figure:: pictures/graph.png
   :scale: 85 %
   :alt: The DAG of graph g
   :align: center

Utility Functions for Graphs (optional)
------------------------------------------
HEMDAG library includes several utility functions to process and analyze graphs as well as I/O functions to import a graph as object of class ``graphNEL`` or to export a graph as object of class ``graphNEL`` in a plain text file (in the classical tupla format). For more details on these functions, please have a look to the `reference manual <https://cran.r-project.org/web/packages/HEMDAG/HEMDAG.pdf>`_.

.. _hem:

Hierarchical Ensemble Methods
================================
First of all, we need to find the root node (i.e. node that is at the top-level of the hierarchy) of the HPO graph ``g``. To do that just type:

.. code-block:: R

    > root <- root.node(g);

in this way we store in the variable ``root`` the root node of the graph ``g``.

Now, we are ready to run any ensemble algorithms implemented in the HEMDAG package.

.. _htd:

HTD-DAG: Hierarchical Top Down for DAG
-----------------------------------------
The HTD-DAG algorithm modifies the flat scores according to the hierarchy of a DAG :math:`G` through a unique run across the nodes of the graph. For a given example :math:`x`, the flat predictions :math:`f(x) = \hat{y}` are hierarchically corrected to :math:`\bar{y}`, by per-level visiting the nodes of the DAG from top to bottom according to the following simple rule:

.. math::

    \bar{y}_i := \left\{
       \begin{array}{lll}
         \hat{y}_i  & {\rm if} \quad i \in root(G) \\
         \min_{j \in par(i)} \bar{y}_j & {\rm if} \quad \min_{j \in par(i)} \bar{y}_j < \hat{y}_i \\
         \hat{y}_i & {\rm otherwise}
       \end{array}
      \right.

The node levels correspond to their maximum path length from the root. To call the HTD-DAG algorithm just type:

.. code-block:: R

    > S.htd <- htd(S.norm, g, root);

Alternatively, we can call the ``htd.vanilla`` function (instead of ``htd``), which it allows to normalize the flat scores matrix ``S`` (according to **maxnorm** or **qnorm** normalization) *on the fly*:

run a normalization method (between **maxnorm** and **qnrom**) *on the fly*:

.. code-block:: R

    > S.htd <- htd.vanilla(S, g, norm=TRUE, norm.type="max.norm");

.. note::

    In ``htd.vanilla``, if ``norm=FALSE`` and ``norm.type=NULL`` the flat scores matrix ``S`` is not normalized.

.. _gpav:

GPAV: Generalized Pool-Adjacent-Violators
--------------------------------------------
Burdakov et al. in [Burdakov06]_ proposed an approximate algorithm, named GPAV, to solve the *isotonic regression* (IR) or *monotonic regression* (MR) problem in its general case (i.e. partial order of the constraints). GPAV algorithm combines both low computational complexity (estimated to be :math:`\mathcal{O}(|V|^2`), where :math:`V` is the number of nodes of the graph) and high accuracy. Formally, given a vector of observed values :math:`\hat{y} \in R^n`, a strictly positive vector of weights :math:`w \in R^n` and a dag :math:`G(V,E)`, GPAV finds the vector of fitted values :math:`\bar{y} \in \mathbb{R}^n` that solves the following convex quadratic program:

.. math::

  \begin{equation}
    \begin{array}{ll}
        \min\limits_{\bar{y}} \quad \sum\limits_{i \in V} w_i (\bar{y}_i - \hat{y}_i)^2 \\
        s.t. \quad \bar{y}_j \geq \bar{y}_i \quad \forall (i,j) \in E
    \end{array}
  \end{equation}

.. [Burdakov06] O. Sysoev, A. Grimvall, and O. Burdakov, Data preordering in generalized pav algorithm for monotonic regression, Journal of Computational Mathematics, vol. 24, no. 6, pp. 771–790, 2006

To call the GPAV algorithm just type:

.. code-block:: R

    > S.gpav <- gpav.over.examples(S.norm, g, W=NULL);

It is worth noting that there is also a parallel version of the GPAV algorithm:

.. code-block:: R

    > S.gpav <- gpav.parallel(S.norm, g, W=NULL, ncores=8);

Similarly to HTD-DAG also for GPAV, we can use the function ``gpav.vanilla`` (instead of ``gpav.over.examples`` or ``gpav.parallel``) to normalize the flat scores matrix ``S`` (according to **maxnorm** or **qnorm** normalization) *on the fly*:

.. code-block:: R

    > S.gpav <- gpav.vanilla(S, g, W=NULL, parallel=FALSE, ncores=8, norm=TRUE, norm.type="maxnorm");

.. _tpr:

TPR-DAG: True Path Rule for DAG
------------------------------------------------
TPR-DAG is a family of algorithms on the basis of the choice of the **bottom-up** step adopted for the selection of *positive* children. Indeed, in their more general form, the TPR-DAG algorithms adopt a two step learning strategy:

    1. in the first step they compute a *per-level bottom-up* visit from leaves to root to propagate *positive* predictions across the hierarchy;
    2. in the second step they compute a *per-level top-down* visit from root to leaves in order to assure the consistency of the predictions. In other word, the :ref:`htd` algorithm is applied.

.. note::

    Levels (both in the first and second step) are defined in terms of the maximum path length from the root node. Please refer to our `BMC Bioinformatics paper <https://doi.org/10.1186/s12859-017-1854-y>`_ for further details.

The *vanilla* TPR-DAG adopts a per-level bottom-up traversal of the DAG to modify the flat predictions :math:`\hat{y}_i` according to the following formula:

.. math::

    \bar{y}_i := \frac{1}{1 + |\phi_i|} (\hat{y}_i + \sum_{j \in \phi_i} \bar{y}_j)

where :math:`\phi_i` are the positive children of :math:`i` (parameter ``positive="children"``).

Different strategies to select the positive children :math:`\phi_i` can be applied:

    1. **threshold-free** strategy (parameter ``bottom="threshold.free"``): the positive nodes are those children that can increment the score of the node :math:`i`, that is those nodes that achieve a score higher than that of their parents:

    .. math::

        \phi_i := \{ j \in child(i) | \bar{y}_j > \hat{y}_i \}

    2. **threshold** strategy (parameter ``bottom="threshold"``): the positive children are selected on the basis of a threshold that can be selected in two different ways:

        a) a unique threshold :math:`\bar{t}` is a priori selected for all nodes to determine the set of positives

        .. math::

            \phi_i := \{ j \in child(i) | \bar{y}_j > \bar{t} \}, \forall i \in V

        For instance if the predictions represent probabilities it could be meaningful to a priori select :math:`\bar{t}=0.5`.

        b) a threshold is selected to maximize some imbalance-aware performance metric :math:`\mathcal{M}` estimated on the training data, as for instance the Fmax or the AUPRC. In other words, the threshold is selected to maximize the measure :math:`\mathcal{M}(j,t)` on the training data for the term :math:`j` with respect to the threshold :math:`t`. The corresponding set of positives for each :math:`i \in V` is:

        .. math::

            \phi_i := \{ j \in child(i) | \bar{y}_j > t_j^*,  t_j^* = \arg \max_{t} \mathcal{M}(j,t) \}

       Internal cross-validation is used to select :math:`t^*_j` within a set of possible thresholds :math:`t \in (0,1)`;

The weighted TPR-DAG version (parameter ``bottom="weighted.threshold.free"``) can be designed by adding a weight :math:`w \in [0,1]` to balance the contribution of the parent node :math:`i` and its positive children :math:`\phi`:

.. math::

    \bar{y}_i := w \hat{y}_i + \frac{(1 - w)}{|\phi_i|} \sum_{j \in \phi_i} \bar{y}_j

If :math:`w=1` no weight is attributed to the children and the TPR-DAG reduces to the HTD-DAG algorithm. If :math:`w=0` only the predictors associated to the children nodes vote to predict node :math:`i`. In the intermediate cases we attribute more importance to the predictor for the node :math:`i` or to its children depending on the values of :math:`w`.

By combining the weighted and the threshold variant, we design the *weighted-threshold* variant (parameter ``bottom="weighted.threshold"``).

All the *vanilla* TPR-DAG variants use the HTD-DAG algorithm in the top-down step (parameter ``topdown="htd"``) to provide ontology-based predictions (i.e. predictions that are coherent with the ontology structure):

.. code-block:: R

    > S.tprTF <- tpr.dag(S.norm, g, root, positive="children", bottomup="threshold.free", topdown="htd");
    > S.tprT  <- tpr.dag(S.norm, g, root, positive="children", bottomup="threshold", topdown="htd", t=0.5);
    > S.tprW  <- tpr.dag(S.norm, g, root, positive="children", bottomup="weighted.threshold.free", topdown="htd", w=0.5);
    > S.tprWT <- tpr.dag(S.norm, g, root, positive="children", bottomup="weighted.threshold", topdown="htd", t=0.5, w=0.5);

DESCENS: Descendants Ensemble Classifier
------------------------------------------------
As shown in [Valentini11]_ for tree-based hierarchies, the contribution of the descendants of a given node decays exponentially with their distance from the node itself and it is straightforward to see that this property also holds for DAG structured taxonomies. To overcame this limitation and in order to enhance the contribution of the most specific nodes to the overall decision of the ensemble we design the ensemble variant DESCENS. The novelty of DESCENS consists in strongly considering the contribution of all the descendants of each node instead of only that of its children (``positive="descendants"``). Therefore DESCENS predictions are more influenced by the information embedded in the leaves nodes, that are the classes containing the most informative and meaningful information from a biological and medical standpoint. DESCENS variants can be designed on the choice of the *positive* descendants :math:`\Delta_i`. The same strategies adopted for the choice of :math:`\phi_i` can be also adopted for the choice of :math:`\Delta_i`, simply by replacing :math:`\phi_i` with :math:`\Delta_i` and :math:`child(i)` with :math:`desc(i)` in the various formulas shown in :ref:`tpr`. Furthermore, we designed a variant specific only for DESCENS, that we named DESCENS-:math:`\tau` (parameter ``bottomup="tau"``). The DESCENS-:math:`\tau` variant balances the contribution between the *positives* children of a node :math:`i` and that of its *positives* descendants excluding its children by adding a weight :math:`\tau \in [0,1]`:

.. math::

    \bar{y}_i := \frac{\tau}{1+|\phi_i|}(\hat{y}_i + \sum_{j \in \phi_i} \bar{y}_j) + \frac{1-\tau}{1+|\delta_i|}(\hat{y}_i + \sum_{j\in \delta_i} \bar{y}_j)

where :math:`\phi_i` are the *positive* children of :math:`i` and :math:`\delta_i=\Delta_i \setminus \phi_i` the descendants of :math:`i` without its children.

If :math:`\tau=1` we consider only the contribution of the *positive* children of :math:`i`; if :math:`\tau=0` only the descendants that are not children contribute to the score, while for intermediate values of :math:`\tau` we can balance the contribution of :math:`\phi_i` and :math:`\delta_i` positive nodes.

.. [Valentini11] G. Valentini, "True Path Rule Hierarchical Ensembles for Genome-Wide Gene Function Prediction," in IEEE/ACM Transactions on Computational Biology and Bioinformatics, vol. 8, no. 3, pp. 832-847, May-June 2011, doi: 10.1109/TCBB.2010.38.

All the DESCENS variants adopt in the second step the HTD-DAG algorithm to assure the consistency of the predictions:

.. code-block:: R

    > S.descensTF  <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="threshold.free", topdown="htd");
    > S.descensT   <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="threshold", topdown="htd", t=0.5);
    > S.descensW   <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="weighted.threshold.free", topdown="htd", w=0.5);
    > S.descensWT  <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="weighted.threshold", topdown="htd", t=0.5, w=05);
    > S.descensTAU <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="tau", topdown="htd", t=0.5);

ISO-TPR: Isotonic Regression for DAG
------------------------------------------------
The ISO-TPR algorithms (parameter ``positive="children"`` and ``topdown="gpav"``) considering the **positive children** in the bottom-up step and adopt GPAV (:ref:`gpav`) instead of HTD-DAG (:ref:`htd`) in the consistency step. The most important feature of the ISO-TPR algorithms is that they maintain the hierarchical constraints by construction by selecting the closest solution (in the least square sense) to the bottom-up predictions that obey the *True Path Rule*:

.. code-block:: R

    > S.isotprTF <- tpr.dag(S.norm, g, root, positive="children", bottomup="threshold.free", topdown="gpav");
    > S.isotprT  <- tpr.dag(S.norm, g, root, positive="children", bottomup="threshold", topdown="gpav", t=0.5);
    > S.isotprW  <- tpr.dag(S.norm, g, root, positive="children", bottomup="weighted.threshold.free", topdown="gpav", w=0.5);
    > S.isotprWT <- tpr.dag(S.norm, g, root, positive="children", bottomup="weighted.threshold", topdown="gpav", t=0.5, w=0.5);

ISO-DESCENS: Isotonic Regression with Descendants Ensemble Classifier
-------------------------------------------------------------------------
The ISO-DESCENS variants (parameter ``positive="descendants"`` and ``topdown="gpav"``) considering the **positive descendants** instead of **positive children** in the bottom-up step and adopt GPAV (instead of the HTD-DAG algorithm) to guarantee the consistency of the predictions:

.. code-block:: R

    > S.isodescensTF  <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="threshold.free", topdown="gpav");
    > S.isodescensT   <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="threshold", topdown="gpav", t=0.5);
    > S.isodescensW   <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="weighted.threshold.free", topdown="gpav", w=0.5);
    > S.isodescensWT  <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="weighted.threshold", topdown="gpav", t=0.5, w=0.5);
    > S.isodescensTAU <- tpr.dag(S.norm, g, root, positive="descendants", bottomup="tau", topdown="gpav", t=0.5);

Obozinski Heuristic Methods
--------------------------------
HEMDAG includes also the three heuristics ensemble methods (And, Max, Or) proposed in [Obozinski08]_:

1. **Max**: reports the largest logistic regression (LR) value of self and all descendants: :math:`p_i = max_{j \in descendants(i)} \hat{p_j}`;

2. **And**: reports the product of LR values of all ancestors and self. This is equivalent to computing the probability that all ancestral terms are "on" assuming that, conditional on the data, all predictions are independent: :math:`p_i = \prod_{j \in ancestors(i)} \hat{p_j}`;

3. **Or**: computes the probability that at least one of the descendant terms is "on" assuming again that, conditional on the data, all predictions are independent: :math:`1 - p_i = \prod_{j \in descendants(i)} (1 - \hat{p_j})`;

.. [Obozinski08] Obozinski G, Lanckriet G, Grant C, M J, Noble WS. Consistent probabilistic output for protein function prediction. Genome Biology. 2008;9:135–142. doi:10.1186/gb-2008-9-s1-s6.

To call Obozinski's heuristic methods, just type:

.. code-block:: R

    > S.max <- obozinski.max(S.norm, g, root);
    > S.and <- obozinski.and(S.norm, g, root);
    > S.or  <- obozinski.or(S.norm, g, root);

Alternatively, the Obozinski's methods can be also called by properly setting the parameter ``heuristic`` of the function ``obozinski.methods``:

.. code-block:: R

    > S.max <- obozinski.methods(S, g, heuristic="max", norm=TRUE, norm.type="maxnorm");
    > S.and <- obozinski.methods(S, g, heuristic="and", norm=TRUE, norm.type="maxnorm");
    > S.or  <- obozinski.methods(S, g, heuristic="or",  norm=TRUE, norm.type="maxnorm");

Hierarchical Constraints Check
==================================
Predictions returned by a flat classifier **do not respect** the *True Path Rule* (since they neglecting the structural information between different ontology terms), whereas the predictions returned by a hierarchical ensemble methods **always obey** the *True Path Rule*. According to this rule a *positive* instance for a class implies *positive* instance for all the ancestors of that class. We can easily check this fact by using the function ``check.hierarchy``. Below (as an example) we check the consistency of the scores corrected according to the HTD-DAG strategy. Of course, all the flat scores corrected with any hierarchical ensemble variants included in HEMDAG, respect the **True Path Rule**. We leave to the reader the freedom to check the consistency of the scores matrix of the remaining 22 hierarchical ensemble variants encompassed in HEMDAG.

.. code-block:: R

    > check.hierarchy(S, g, root)$status
    [1] "NOTOK"

    > check.hierarchy(S.htd, g, root)$status
    [1] "OK"

Performance Evaluation
==========================
To know the behavior of the hierarchical ensemble methods, the HEMDAG library provides both *term-centric* and *protein-centric* performance metrics:

- ``AUPRC``: area under the precision-recall curve;
- ``AUROC``: area under the ROC curve;
- ``Fmax`` : maximum hierarchical F-score [Jiang2016]_;
- ``PXR``  : precision at different recall levels;

.. note::
    a) HEMDAG allows to compute all the aforementioned performance metrics either **one-shot** or **averaged** across k fold. Depending on the dataset size, the metrics ``Fmax`` and ``PXR`` could take a while to finish. Please refer to HEMDAG `reference manual <https://cran.r-project.org/web/packages/HEMDAG/HEMDAG.pdf>`_  for further information about the input arguments of these functions.
    b) For computing the *term-centric* metrics (``AUROC``, ``AUPRC`` and ``PXR``), HEMDAG makes use of the R package *precrec* (`link <https://CRAN.R-project.org/package=precrec>`__).

.. [Jiang2016] Y. Jiang et al., An expanded evaluation of protein function prediction methods shows an improvement in accuracy, Genome Biology, vol. 17, p. 184, 2016

Load the Annotation Matrix
------------------------------
To compare the hierarchical ensemble methods against the flat approach, we need to load the annotation matrix:

.. code-block:: R

    > data(labels);

With the above command we loaded the annotations table ``L``, that is a named ``100 X 23`` matrix. Rows correspond to genes (``Entrez GeneID``) and columns to HPO terms/classes. ``L[i, j] = 1`` means that the gene ``i`` belong to class ``j``, ``L[i, j] = 0`` means that the gene ``i`` does not belong to class ``j``.

Flat vs Hierarchical
------------------------
Before computing performance metrics we should remove the root node from the annotation matrix, the flat scores matrix and the hierarchical scores matrix. Indeed, it does not make sense to take into account the predictions of the root node, since it is a *fake* node added to the ontology for practical reasons (e.g. some graph-based software may require a single root node to work). In R this can be accomplished in one line of code.

.. code-block:: R

    ## remove root node from annotation matrix
    > if(root %in% colnames(L))
    +    L <- L[,-which(colnames(L)==root)];

    ## remove root node from the normalized flat scores matrix
    > if(root %in% colnames(S.norm))
    +    S.norm <- S.norm[,-which(colnames(S.norm)==root)];

    ## remove root node from hierarchical scores matrix (eg S.htd)
    > if(root %in% colnames(S.htd))
    +    S.htd <- S.htd[,-which(colnames(S.htd)==root)];

Now we can compare the flat approach RANKS versus the HTD-DAG strategy, by averaging (for instance) the performance across 3 folds:

.. code-block:: R

    ## RANKS
    > prc.flat  <- auprc.single.over.classes(L, S.norm, folds=3, seed=23);
    > auc.flat  <- auroc.single.over.classes(L, S.norm, folds=3, seed=23);
    > pxr.flat  <- precision.at.given.recall.levels.over.classes(L, S.norm, recall.levels=seq(from=0.1, to=1, by=0.1), folds=3, seed=23);
    > fmax.flat <- compute.fmax(L, S.norm, n.round=3, verbose=FALSE, b.per.example=TRUE, folds=3, seed=23);

    ## HTD-DAG
    > prc.htd  <- auprc.single.over.classes(L, S.htd, folds=3, seed=23);
    > auc.htd  <- auroc.single.over.classes(L, S.htd, folds=3, seed=23);
    > pxr.htd  <- precision.at.given.recall.levels.over.classes(L, S.htd, recall.levels=seq(from=0.1, to=1, by=0.1), folds=3, seed=23);
    > fmax.htd <- compute.fmax(L, S.htd, n.round=3, verbose=FALSE, b.per.example=TRUE, folds=3, seed=23);

By looking at the results, it easy to see that the HTD-DAG outperforms the flat classifier RANKS:

.. code-block:: R

   ## AUC performance: RANKS VS HTD-DAG
    > auc.flat$average
    [1] 0.8297
    > auc.htd$average
    [1] 0.8336

    ## PRC performance: RANKS VS HTD-DAG
    > prc.flat$average
    [1] 0.4333
    > prc.htd$average
    [1] 0.4627

    ## Fmax performance: RANKS VS HTD-DAG
    > fmax.flat$average
        P      R      S      F    avF      A      T
    0.5042 0.8639 0.4485 0.6368 0.5269 0.6612 0.5720
    > fmax.htd$average
        P      R      S      F    avF      A      T
    0.5576 0.7745 0.6519 0.6484 0.5617 0.7521 0.6487

    ## PXR: RANKS VS HTD-DAG
    > pxr.flat$average
       0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1
    0.5821 0.5821 0.5821 0.5531 0.5531 0.4483 0.4388 0.4388 0.4388 0.4388
    > pxr.htd$average
       0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1
    0.6218 0.6218 0.6218 0.5941 0.5941 0.4798 0.4668 0.4668 0.4668 0.4668

.. note::
    HTD-DAG is the simplest ensemble approach among those available. HTD-DAG strategy makes flat scores consistent with the hierarchy by propagating from top to bottom the negative predictions. Hence, in the worst case might happen that the predictions at leaves nodes are all negatives. Other ensemble algorithms, such as GPAV and TPR-DAG (and variants) should lead to better improvements.

Tuning of Hyper-Parameter(s)
===============================
14 out of 18 of the TPR-DAG hierarchical algorithms are parametric. Instead of use a priori selected threshold (as done in :ref:`tpr` and variants), we can tune the hyper-parameter(s) of the parametric variants through the function ``tpr.dag.cv``. The hyper-parameter(s) can be maximize on the basis of ``AUPRC`` (parameter ``metric="prc"``) or ``Fmax`` (parameter ``metric="fmax"``). Below, as an example, we maximize the threshold of the parametric variant ISO-TPR-Threshold (``isotprT``) on the basis of ``AUPRC`` metric.

.. code-block:: R

    > threshold <- seq(0.1, 0.9, 0.1);

    > S.isotprT <- tpr.dag.cv(S, g, ann=L, norm=TRUE, norm.type="maxnorm", positive="children",
                              bottomup="threshold", topdown="gpav", W=NULL, parallel=FALSE,
                              ncores=1, threshold=threshold, weight=0, kk=3, seed=23,
                              metric="prc", n.round=NULL);

    ## stdout
    maxnorm normalization: done
    training fold:  1   top prc avg found:  0.4536119   best threshold: 0.1
    training fold:  1   top prc avg found:  0.4592147   best threshold: 0.4
    training fold:  2   top prc avg found:  0.2190192   best threshold: 0.1
    training fold:  2   top prc avg found:  0.2193331   best threshold: 0.6
    training fold:  2   top prc avg found:  0.2208776   best threshold: 0.7
    training fold:  3   top prc avg found:  0.8148121   best threshold: 0.1
    tpr-dag correction done

Evaluating ``isotprT`` by computing *term-* and *protein-* centric performance (always averaging the performance across 3 folds), it easy to see how this ensemble variant outperform both the flat classifier RANKS and the hierarchical algorithm HTD-DAG:

.. code-block:: R

    ## remove root node before computing performance
    > if(root %in% colnames(S.isotprT))
    +    S.isotprT <- S.isotprT[,-which(colnames(S.isotprT)==root)];

    > prc.isotprT  <- auprc.single.over.classes(L, S.isotprT, folds=3, seed=23);
    > auc.isotprT  <- auroc.single.over.classes(L, S.isotprT, folds=3, seed=23);
    > pxr.isotprT  <- precision.at.given.recall.levels.over.classes(L, S.isotprT, recall.levels=seq(from=0.1, to=1, by=0.1), folds=3, seed=23);
    > fmax.isotprT <- compute.fmax(L, S.isotprT, n.round=3, verbose=FALSE, b.per.example=TRUE, folds=3, seed=23);

    ## AUC performance: RANKS VS HTD-DAG vs isotprT
    > auc.flat$average
    [1] 0.8297
    > auc.htd$average
    [1] 0.8336
    > auc.isotprT$average
    [1] 0.8446

    ## PRC performance: RANKS VS HTD-DAG vs isotprT
    > prc.flat$average
    [1] 0.4333
    > prc.htd$average
    [1] 0.4627
    > prc.isotprT$average
    [1] 0.5346

    ## Fmax performance: RANKS VS HTD-DAG vs isotprT
    > fmax.flat$average
        P      R      S      F    avF      A      T
    0.5042 0.8639 0.4485 0.6368 0.5269 0.6612 0.5720
    > fmax.htd$average
        P      R      S      F    avF      A      T
    0.5576 0.7745 0.6519 0.6484 0.5617 0.7521 0.6487
    > fmax.isotprT$average
        P      R      S      F    avF      A      T
    0.5896 0.8306 0.5283 0.6896 0.6106 0.7066 0.6340

    ## PXR: RANKS VS HTD-DAG vs isotprT
    > pxr.flat$average
       0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1
    0.5821 0.5821 0.5821 0.5531 0.5531 0.4483 0.4388 0.4388 0.4388 0.4388
    > pxr.htd$average
       0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1
    0.6218 0.6218 0.6218 0.5941 0.5941 0.4798 0.4668 0.4668 0.4668 0.4668
    > pxr.isotprT$average
       0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1
    0.6848 0.6848 0.6848 0.6697 0.6697 0.5417 0.5027 0.5027 0.5027 0.5027

By properly setting the parameters ``positive``, ``bottomup`` and ``topdown`` of the function ``tpr.dag.cv``, it is easy to make experiments with all the 18 TPR-DAG ensemble variants. For further details on the other input arguments of the function ``tpr.dag.cv``, please refer to the `reference manual <https://cran.r-project.org/web/packages/HEMDAG/HEMDAG.pdf>`_.

.. note::

    Note that tuning the hyper-parameter(s) of the ensemble variants on the basis of ``Fmax`` might involve high running time (due to the nature itself of the ``Fmax`` metric).

Hold-out Functions
===================
For all the hierarchical ensemble algorithms encompassed in the HEMDAG library there is also a corresponding hold-out version. The hold-out functions respect to the *vanilla* ones, require in input a vector of integer numbers corresponding to the indexes of the elements (rows) of the scores matrix ``S`` to be used in the test set (parameter ``testIndex``). The hold-out ensemble functions included in HEMDAG are:

    * ``htd.holdout``;
    * ``gpav.holdout``;
    * ``tpr.dag.holdout``;
    * ``obozinski.holdout``;

For the sake of space we do not show here experiments by using the hold-out version of the hierarchical functions. Please refer to the `reference manual <https://cran.r-project.org/web/packages/HEMDAG/HEMDAG.pdf>`_, for further details on these functions.
