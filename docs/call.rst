.. role:: raw-html(raw)
   :format: html

.. _call:

====================
HEMDAG Playground
====================
Here we explain how to play with the ensemble algorithms of the HEMDAG family in both cross-validated and hold-out experiments.

HEMDAG can in principle boost the predictions of any flat learning methods by reconciling the flat predictions with the topology of the underlying ontology. Then, before using HEMDAG we need the following *ingredients*:

1) the labels matrix ``M`` representing the proteins' annotations to functional terms;
2) the graph ``g`` representing the hierarchy of the functional terms;
3) the flat scores matrix ``S`` representing a score or a probability that a gene/protein belonging to a given functional term;

For instance, to build dataset for the GO term prediction you can use the following `pipeline <https://github.com/marconotaro/godata-pipe>`__; instead to obtain flat predictions you can use the `shogun library <https://www.shogun-toolbox.org/>`__ or the `caret package <https://topepo.github.io/caret/>`__ or any other software returning a score or a probability that a protein belong to a functional term.

.. note::

    HEMDAG builds upon flat predictions. Consequently, the predictions returned by the adopted learning method **must violate** the hierarchical relationships between ontology terms, otherwise the application of HEMDAG is meaningless.


Call HEMDAG for time-lapse hold-out experiments
===================================================
For the sake of the simplicity, here we use pre-built datasets for the organism *Drosophila melanogaster* (DROME) on the GO domain molecular function (MF). The datasets were built by using a classical time-lapse hold-out procedure. More precisely, we used the annotations of an old GO release (December 2017) to predict the protein functions of a more recent GO release (June 2020). The graph and the annotation matrix was built by adopting the following `pipeline <https://github.com/marconotaro/godata-pipe#build-dataset-for-a-time-split-hold-out-procedure>`__. The flat scores matrix was obtained by using the R interface of the machine learning library `LiblineaR <https://CRAN.R-project.org/package=LiblineaR>`__ with the default parameter settings.


Download data
----------------
All the datasets can be downloaded by using the following commands:

.. code-block:: bash

    mkdir -p ~/hemdag-playground/data/
    cd ~/hemdag-playground/data/
    curl -Ss https://github.com/marconotaro/hemdag/tree/master/docs/playground/data/ho |  grep -oP '(?<=href=").*?(?=\">)' | grep '.rda$' | perl -pe 's/blob\///' | perl -pe 's/^/https\:\/\/raw.githubusercontent.com/' | wget -nc -i -

With command above we download the following datasets:

``7227.drome.go.mf.ann.20dec17-16jun20.rda``: the annotation matrix;
``7227.drome.go.mf.dag.20dec17-16jun20.rda``: the graph;
``7227.drome.go.mf.testindex.20dec17-16jun20.rda``: the indexes of the examples of the test set;
``7227.drome.go.mf.scores.svmlinear.holdout.rda``: the flat scores matrix;


HEMDAG hold-out script
-------------------------
To call an HEMDAG algorithm on a time-lapse hold-out dataset you can execute the script shown below. To download the script:

.. code-block:: bash

    mkdir -p ~/hemdag-playground/script/
    cd ~/hemdag-playground/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/hemdag-call.R

Before executing the script be sure to have correctly installed the latest version of the HEMDAG package (and all its dependencies -- see :ref:`installation`) and the package `optparse <https://cran.r-project.org/web/packages/optparse/>`__.

.. note::

    the code lines highlighted in yellow are the lines that differ between the script to call HEMDAG in hold-out experiments and the script to call HEMDAG in cross-validated experiments. See :ref:`hemdagcv`.

.. literalinclude:: playground/script/hemdag-call.R
    :language: R
    :linenos:
    :emphasize-lines: 111,112,120,147-159


Arguments Explanation
~~~~~~~~~~~~~~~~~~~~~~~~~~~
For the usage of the script, type in the bash under the ``~/hemdag-playground/script/`` folder:

.. code-block:: bash

    Rscript hemdag-call.R -h

For a detailed description of the input arguments of the ``hemdag-call.R`` script, please refer to the input variable of the functions ``(gpav|htd|tpr.dag).holdout`` in the `HEMDAG reference manual <https://cran.r-project.org/web/packages/HEMDAG/HEMDAG.pdf>`__.

Parametric-free variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The arguments to be set to call the parametric-free HEMDAG algorithms are:

* ``-b (--bottomup)`` :raw-html:`&rarr;` ``none``
* ``-t (--topdown)``  :raw-html:`&rarr;` ``gpav|htd``

Optional: GPAV can also be run in parallel simply by setting the flag ``-l (--parallel)`` and the number of cores ``-n (--cores)``.

Parametric variants
~~~~~~~~~~~~~~~~~~~~~~~~
The arguments to be set to call the parametric HEMDAG algorithms are:

* ``-b (--bottomup)`` :raw-html:`&rarr;` ``threshold.free|threshold|weighted.threshold.free|weighted.threshold|tau``
* ``-t (--topdown)``  :raw-html:`&rarr;` ``gpav|htd``
* ``-c (--cut-off)``  :raw-html:`&rarr;` ``0.5|"seq(from=0.1, to=0.9, by=0.1)"``. Note the use of double quotes for range of thresholds
* ``-w (--weight)``   :raw-html:`&rarr;` ``0.5|"seq(from=0.1, to=0.9, by=0.1)"``. Note the use of double quotes for range of thresholds

If a range of thresholds (or weights) is selected, the hyper-parameters are tuned on the basis of imbalance-aware performance metrics estimated on the training data -- ``-m (--metric)`` :raw-html:`&rarr;` ``auprc|fmax``. By default the number of folds ``-k (--fold)`` is set to 5 and the seed ``-s (--seed)`` is set to 23. In addition, if ``-m (--metric)`` :raw-html:`&rarr;` ``fmax`` the parameter ``-r (--round)`` can be used to select the number of rounding digits to be applied for choosing the best Fmax (by default is set to 3).

Optional arguments
~~~~~~~~~~~~~~~~~~~~~~~~~
In case the flat learning method returns a score and not a probability, the flat scores matrix must be normalized before running HEMDAG. On the contrary, if the flat classifier already returns a probability there is no needed to normalize the flat scores matrix, since the flat scores can be directly compared with the hierarchical ones. To normalize the flat scores matrix we must "activate" the flag ``-z (--norm)`` (by default the flag ``-z`` is deactivate) and we need to choose a type of normalization (``-y (--normtype)`` :raw-html:`&rarr;` ``maxnorm|qnorm``).

.. note::

    The name of the chosen normalization is saved in the *rda* output file name.


HEMDAG Call
---------------
Calling different HEMDAG ensemble algorithms it is straightforward as shown in the examples below. The hierarchical scores matrices are stored in the folder ``../res/`` with the name of the chosen HEMDAG algorithm in the output *.rda* file name. The HEMDAG elapsed time is printed on the shell.

.. note::

    #. If the ``hemdag-call.R`` script is called without specified any input arguments, by default the ``isodescensTAU`` algorithm with the tuning of parameter on the basis of AUPRC metric is executed.

    #. The tuning of the hyper-parameters can take from few minutes up to few hours depending on the size of the dataset and on the adopted evaluation metric (Fmax is slower than AUPRC). In the examples shown below the tuning of the hyper-parameter takes few minutes.

``GPAV`` call (parallel version):

.. code-block:: bash

    Rscript hemdag-call.R -o 7227.drome -d mf -b none -t gpav -l -n 12


``isotprTF`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227.drome -d mf -p children -b threshold.free -t gpav -l -n 12


``isotprW`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227.drome -d mf -p children -b weighted.threshold.free -t gpav -w "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12


``isodescensTF`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227.drome -d mf -p descendants -b threshold.free -t gpav -l -n 12


``isodescensW`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227.drome -d mf -p descendants -b weighted.threshold.free -t gpav -w "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12


``isodescensTAU`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227.drome -d mf -p descendants -b tau -t gpav -c "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12


Check Hierarchical Constraints
--------------------------------
All the HEMDAG scores matrices respect the hierarchical constraints imposed by the underlying ontology. As an example, the script below checks that all the 6 HEMDAG matrices obtained with the command above do not violate the between-term relationships in the GO hierarchy. Refer to :ref:`conscheck` for  more details.

.. literalinclude:: playground/script/hemdag-constraints-check.R
    :language: R
    :linenos:
    :emphasize-lines: 7,6,20,21


To download the script:

.. code-block:: bash

    mkdir -p ~/hemdag-playground/script/
    cd ~/hemdag-playground/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/hemdag-constraints-check.R


To call the script:

.. code-block:: R

    Rscript hemdag-constraints-check.R

    ## stdout
    drome mf svmlinear+gpav check passed :)
    drome mf svmlinear+isotprTF check passed :)
    drome mf svmlinear+isotprW check passed :)
    drome mf svmlinear+isodescensTF check passed :)
    drome mf svmlinear+isodescensW check passed :)
    drome mf svmlinear+isodescensTAU check passed :)


Evaluate HEMDAG
-------------------
To evaluate the generalization performance of the hierarchical ensemble methods, you must load in the R environment the annotation matrix and the hierarchical scores matrix. Then you can use the *term-centric* and/or *protein-centric* evaluation metrics provided by the HEMDAG package. For further details refer to :ref:`eval`.

.. literalinclude:: playground/script/hemdag-perf-eval.R
    :language: R
    :linenos:

To download the script:

.. code-block:: bash

    mkdir -p ~/hemdag-playground/script/
    cd ~/hemdag-playground/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/hemdag-perf-eval.R


To call the script on single dataset:

.. code-block:: bash

    Rscript hemdag-perf-eval.R -o 7227.drome -d mf -f svmlinear -a gpav


To evaluate more datasets in parallel, you can use the perl script below to generate multiple evaluation call:

.. literalinclude:: playground/script/call-hemdag-perf-eval.pl
    :linenos:
    :language: perl

then you can save them in a bash file and execute it:

.. code-block:: bash

    ## download perl script
    mkdir -p ~/hemdag-playground/script/
    cd ~/hemdag-playground/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/call-hemdag-perf-eval.pl

    ## generate multiple evaluation calls
    perl call-hemdag-perf-eval.pl > hemdag-perf-eval.sh

    ## run multiple evaluations in parallel
    bash hemdag-perf-eval.sh > out


.. _hemdagcv:

Call HEMDAG for cross-validated experiments
===============================================

what to do here:

    * highlight the differences respect to holdout procedure
    * load data set (drome mf ranger)
    * add a flag in the above scripts to call HEMDAG on cv/ho experiments :raw-html:`&rarr;` remove yellow code lines

