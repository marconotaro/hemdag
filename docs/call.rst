.. role:: raw-html(raw)
   :format: html

.. _call:

===========================================
HEMDAG programmatic call and evaluation
===========================================
Here we explain how to apply the ensemble algorithms of the HEMDAG family in both cross-validated and hold-out experiments.

HEMDAG can in principle boost the predictions of any flat learning method by reconciling the flat predictions with the topology of the underlying ontology. Hence, to run HEMDAG we need the following *ingredients*:

1) the label matrix ``M`` representing the protein annotations to functional terms;
2) the graph ``g`` representing the hierarchy of the functional terms;
3) the flat score matrix ``S`` representing a score or a probability that a gene/protein belongs to a given functional term;

To build the graph, the label matrix and the protein-protein interaction network you can use this `pipeline <https://github.com/marconotaro/godata-pipe>`__.
Instead, to obtain the flat score matrix you can use the `shogun library <https://www.shogun-toolbox.org/>`__ or the `caret package <https://topepo.github.io/caret/>`__ or any other software returning a score or a probability that a protein belongs to a functional term.

.. note::

    HEMDAG is built upon flat predictions. HEMDAG corrects all the violations of the hierarchical relationships between ontology terms.

.. _hemdagscript:

HEMDAG Call Script
======================
To call any hierarchical ensemble algorithm of the HEMDAG family on either a time-lapse hold-out or a cross-validated dataset you can execute the following script:

.. literalinclude:: playground/script/hemdag-call.R
    :language: R
    :linenos:

You can download the script as follow:

.. code-block:: bash

    mkdir -p ~/hemdag/script/
    cd ~/hemdag/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/hemdag-call.R

Before executing the script be sure to have correctly installed the latest version of the HEMDAG package (and all its dependencies -- see :ref:`installation`) and the package `optparse <https://cran.r-project.org/web/packages/optparse/>`__.


.. note::

    #. The output hierarchical score matrix of the called HEMDAG algorithm (whose name is saved in the output *.rda* file name) is stored in the folder ``~/hemdag/res/(ho|cv)`` depending on whether you chose to execute HEMDAG on either hold-out (ho) or cross-validated (cv) datasets. The HEMDAG elapsed time is printed on the shell.
    #. By default, if no inputs parameters are specified in ``hemdag-call.R``, the script executes the ``isodescensTAU`` algorithm on the hold-out dataset by tuning the parameter ``tau`` on the basis of AUPRC.
    #. The tuning of the hyper-parameters can take from few minutes up to few hours depending on the size of the dataset and on the adopted evaluation metric (Fmax is slower than AUPRC).

Arguments Explanation
-------------------------
For the usage of the script, type in the shell under the ``~/hemdag/script/`` folder:

.. code-block:: bash

    Rscript hemdag-call.R -h


For a detailed description of the input arguments *positive, bottomup, topdown, threshold, weight, metric, round, seed, fold, parallel, cores, norm, normtype*, please refer to the description of the input variables of the functions ``(gpav|htd|tpr.dag).(holdout|cv)`` in the `HEMDAG reference manual <https://cran.r-project.org/web/packages/HEMDAG/HEMDAG.pdf>`__.

Parametric-free arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~
To call a parametric-free HEMDAG algorithm the main required arguments are:

* ``-b (--bottomup)`` :raw-html:`&rarr;` ``none``
* ``-t (--topdown)``  :raw-html:`&rarr;` ``gpav|htd``

.. note::

    GPAV can also be run in parallel simply by using the flag ``-l (--parallel)`` and by setting the number of cores ``-n (--cores)``.

Parametric arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~
To call a parametric HEMDAG algorithm the main required arguments are:

* ``-b (--bottomup)`` :raw-html:`&rarr;` ``threshold.free|threshold|weighted.threshold.free|weighted.threshold|tau``
* ``-t (--topdown)``  :raw-html:`&rarr;` ``gpav|htd``
* ``-c (--cut-off)``  :raw-html:`&rarr;` ``0.5|"seq(from=0.1, to=0.9, by=0.1)"``. Note the use of double quotes for the range of thresholds
* ``-w (--weight)``   :raw-html:`&rarr;` ``0.5|"seq(from=0.1, to=0.9, by=0.1)"``. Note the use of double quotes for the range of thresholds

If a range of thresholds (or weights) is selected, the hyper-parameters are tuned on the basis of imbalance-aware performance metrics estimated on the training data -- ``-m (--metric)`` :raw-html:`&rarr;` ``auprc|fmax``. By default the number of folds ``-k (--fold)`` is set to 5 and the seed ``-s (--seed)`` for the random generator is set to 23. Furthermore, if ``-m (--metric)`` :raw-html:`&rarr;` ``fmax`` the parameter ``-r (--round)`` can be used to select the number of rounding digits to be applied for choosing the best Fmax (by default is set to 3).

Additional arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following arguments are dataset-specific:

* ``-o (--organism)`` specifies the organism name (in the form <taxon>_<name>);
* ``-d (--domain)``   specifies the GO domain: bp (biological process), mf (molecular function), cc(cellular component);
* ``-e (--exptype)``  specifies the type of dataset where running HEMDAG: ho (holdout) or cv (cross-validated);
* ``-f (--flat)``     specifies the name of the flat classifier. In case the flat learning method returns a score and not a probability, the flat score matrix must be normalized before running HEMDAG. On the contrary, if the flat classifier already returns a probability there is no needed to normalize the flat score matrix, since the flat scores can be directly compared with the hierarchical ones. To normalize the flat score matrix we must "activate" the flag ``-z (--norm)`` (by default the flag ``-z`` is deactivate) and we need to choose a type of normalization (``-y (--normtype)`` :raw-html:`&rarr;` ``maxnorm|qnorm``). The name of the chosen normalization is stored in the *rda* output file name.


.. _timesplit:

Time-lapse hold-out experiments
====================================
Here, to show how to use HEMDAG in time-lapse hold-out experiments, we use a pre-built dataset of the organism *Drosophila melanogaster* (DROME) and, for simplicity, we consider the annotations of the GO domain molecular function (MF). To build the dataset we used the annotations of an old GO release (December 2017) as training set and the annotations of a more recent GO release (June 2020) as test set. The graph and the annotation matrix was built by adopting the `pipeline <https://github.com/marconotaro/godata-pipe#build-dataset-for-a-time-split-hold-out-procedure>`__. The flat score matrix was obtained by using the R interface of the machine learning library `LiblineaR <https://CRAN.R-project.org/package=LiblineaR>`__ with the default parameter settings. For further details on the dataset, please refer to *HEMDAG: a family of modular and scalable hierarchical ensemble methods to improve Gene Ontology term prediction (submitted to Bioinformatics)*.

Download Data
----------------
All the required *.rda* files can be downloaded by using the following commands, by exploiting the beauty and power of the non-greedy positive lookbehind regex :raw-html:`&#128521;`:

.. code-block:: bash

    mkdir -p ~/hemdag/data/ho/
    cd ~/hemdag/data/ho/
    curl -Ss https://github.com/marconotaro/hemdag/tree/master/docs/playground/data/ho |  grep -oP '(?<=href=").*?(?=\">)' | grep '.rda$' | perl -pe 's/blob\///' | perl -pe 's/^/https\:\/\/raw.githubusercontent.com/' | wget -nc -i -

With the command above, we download the following datasets:

* ``7227_drome_go_mf_ann_20dec17_16jun20.rda``: the annotation matrix;
* ``7227_drome_go_mf_dag_20dec17_16jun20.rda``: the graph;
* ``7227_drome_go_mf_testindex_20dec17_16jun20.rda``: the indexes of the examples of the test set;
* ``7227_drome_go_mf_scores_svmlinear_holdout.rda``: the flat score matrix;

.. _hemdagcall:

Programmatic Call
-------------------
Below we show some examples of how to call HEMDAG in time-lapse hold-out experiments, but we leave the user the freedom to experiment with any another ensemble algorithm of the HEMDAG family.

.. note::

    #. the ``hemdag-call.R`` script must be called in ``~/hemdag/script/``;
    #. for the examples shown below, the tuning of hyper-parameters takes few minutes;
    #. the output HEMDAG score matrices are stored in ``~/hemdag/res/ho/``.


``GPAV`` call (parallel version):

.. code-block:: bash

    Rscript hemdag-call.R -o 7227_drome -d mf -e ho -f svmlinear -b none -t gpav -l -n 12


``isotprTF`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227_drome -d mf -e ho -f svmlinear  -p children -b threshold.free -t gpav -l -n 12


``isotprW`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227_drome -d mf -e ho -f svmlinear  -p children -b weighted.threshold.free -t gpav -w "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12


``isodescensTF`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227_drome -d mf -e ho -f svmlinear  -p descendants -b threshold.free -t gpav -l -n 12


``isodescensW`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227_drome -d mf -e ho -f svmlinear -p descendants -b weighted.threshold.free -t gpav -w "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12


``isodescensTAU`` call:

.. code-block:: bash

    Rscript hemdag-call.R -o 7227_drome -d mf -e ho -f svmlinear -p descendants -b tau -t gpav -c "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12

.. _hemdagcheck:

Check Hierarchical Constraints
--------------------------------
All the HEMDAG score matrices respect the hierarchical constraints imposed by the underlying ontology. The script below checks that all the 6 HEMDAG matrices obtained with the commands shown above, do not violate the between-term relationships in the GO MF hierarchy. For further details refer to :ref:`conscheck`.


.. literalinclude:: playground/script/hemdag-checker.R
    :language: R
    :linenos:


To download and use the performance evaluation script:

.. code-block:: bash

    ## download
    mkdir -p ~/hemdag/script/
    cd ~/hemdag/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/hemdag-checker.R

    ## call
    Rscript hemdag-checker.R -e ho

    ## example of stdout
    drome mf svmlinear+gpav check passed :)
    drome mf svmlinear+isotprTF check passed :)
    drome mf svmlinear+isotprW check passed :)
    drome mf svmlinear+isodescensTF check passed :)
    drome mf svmlinear+isodescensW check passed :)
    drome mf svmlinear+isodescensTAU check passed :)


You can customize the R script ``hemdag-checker.R`` by extending the vectors *orgs, flats, algs, onts* with the values of your interest.

.. _hemdageval:

Evaluation
-------------------
To evaluate the generalization performance of HEDMAG in time-lapse hold-out experiments, you can use the *term-centric* and/or *protein-centric* evaluation metrics provided by the HEMDAG package itself. For further details on the implemented performance metric refer to section :ref:`eval` of the HEMDAG tutorial.

.. literalinclude:: playground/script/hemdag-eval.R
    :language: R
    :linenos:

To download and use the HEMDAG evaluation script:

.. code-block:: bash

    ## download
    mkdir -p ~/hemdag/script/
    cd ~/hemdag/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/hemdag-eval.R

    ## call
    Rscript hemdag-eval.R -o 7227_drome -d mf -e ho -f svmlinear -a gpav


Chunk evaluation (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The above R call evaluates the performance of an HEMDAG algorithm just on a single dataset. Since HEMDAG can be virtually applied on top of any flat classifier, the Perl script below generates "chunks" of HEMDAG evaluation calls.

.. note::

    The parameter regulating the number of chunk is ``$m``, set by default to 12. Increase (resp. decrease) this value to rise (reduce) the number of HEMDAG evaluation calls to be executed in parallel.

.. note: use perl6 and not perl (in language) to highlight syntax, because the perl highlighting fails due to the % modulo operator (probably the module operator causes a parsing error pygments)
.. literalinclude:: playground/script/hemdag-chunk.pl
    :language: perl6
    :linenos:

To download and use the Perl script that generates "chunks" of HEMDAG evaluation calls:

.. code-block:: bash

    ## download the perl script
    mkdir -p ~/hemdag/script/
    cd ~/hemdag/script/
    wget https://raw.githubusercontent.com/marconotaro/hemdag/master/docs/playground/script/hemdag-chunk.pl

    ## generate chunk evaluation calls
    perl hemdag-chunk.pl -e ho -c 6 > hemdag-ho-eval.sh

    ## evaluate HEMDAG in chunks
    bash hemdag-ho-eval.sh > out


You can customize the Perl script ``hemdag-chunk.pl`` by extending the arrays *@orgs, @flats, @algs, @onts* with the values of your interest. For instance, by setting ``my @onts=  qw(bp mf cc)``, the call ``perl call-hemdag-perf-eval.pl -e ho -c 12`` returns in output 2 chunks of evaluation calls, the first made of 12 calls and the second one of 6 calls. Modify the script to see the output printed on the shell :raw-html:`&#128515;`.


Cross-validated experiments
===============================================
Here, to show how to use HEMDAG in cross-validated experiments, we use a pre-built dataset of the organism *Drosophila melanogaster* (DROME) that covers the annotations of the GO domain molecular function (MF). The graph and protein-GO term associations belong to the GO release of December 2017. The graph and the annotation matrix was built by adopting the following `pipeline <https://github.com/marconotaro/godata-pipe##build-dataset-for-a-cross-validation-procedure>`__. The flat score matrix was obtained by using the random forest as flat learning method (model *ranger* in the R library `caret package <https://topepo.github.io/caret/>`__ with the default parameter settings). For further details on the dataset, please refer to *HEMDAG: a family of modular and scalable hierarchical ensemble methods to improve Gene Ontology term prediction (submitted to Bioinformatics)*.

Download Data
-----------------
All the required *.rda* files can be downloaded by using the following commands:

.. code-block:: bash

    mkdir -p ~/hemdag/data/cv/
    cd ~/hemdag/data/cv/
    curl -Ss https://github.com/marconotaro/hemdag/tree/master/docs/playground/data/cv |  grep -oP '(?<=href=").*?(?=\">)' | grep '.rda$' | perl -pe 's/blob\///' | perl -pe 's/^/https\:\/\/raw.githubusercontent.com/' | wget -nc -i -

.. note::

    Note the change of the last directory from ``ho/`` to ``cv/`` compared to the data downloaded in the :ref:`timesplit`.

With the command above, we download the following datasets:

* ``7227_drome_go_mf_ann_20dec17.rda``: the annotation matrix;
* ``7227_drome_go_mf_dag_20dec17.rda``: the graph;
* ``7227_drome_go_mf_scores_pearson_100_feature_ranger_5fcv.rda``: the flat score matrix;


Programmatic Call
---------------------
To execute any HEMDAG algorithm on cross-validated datasets, you must simply replace ``-e ho`` with ``-e cv`` in the various calls shown in section :ref:`hemdagcall` for the time-split hold-out experiments.

.. note::

    #. the ``hemdag-call.R`` script must be called in ``~/hemdag/script/``;
    #. for the examples shown below the tuning of hyper-parameters takes around one hour;
    #. note that the flat classifier used here is not the *svm* (``-f svmlinear``), but the *random forest* (``-f ranger``);
    #. the output HEMDAG score matrices are stored in ``~/hemdag/res/cv/`` (note the shift of the last directory);

For instance, to call the 6 HEMDAG on cross-validated datasets just type:

.. code-block:: bash

    ## GPAV
    Rscript hemdag-call.R -o 7227_drome -d mf -e cv -f ranger -b none -t gpav -l -n 12

    ## isotprTF
    Rscript hemdag-call.R -o 7227_drome -d mf -e cv -f ranger -p children -b threshold.free -t gpav -l -n 12

    ## isotprW
    Rscript hemdag-call.R -o 7227_drome -d mf -e cv -f ranger -p children -b weighted.threshold.free -t gpav -w "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12

    ## isodescensTF
    Rscript hemdag-call.R -o 7227_drome -d mf -e cv -f ranger -p descendants -b threshold.free -t gpav -l -n 12

    ## isodescensW
    Rscript hemdag-call.R -o 7227_drome -d mf -e cv -f ranger -p descendants -b weighted.threshold.free -t gpav -w "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12

    ## isodescensTAU
    Rscript hemdag-call.R -o 7227_drome -d mf -e cv -f ranger -p descendants -b tau -t gpav -c "seq(from=0.1, to=0.9, by=0.1)" -m auprc -s 23 -k 5 -l -n 12


Check Hierarchical Constraints
---------------------------------
To check that HEMDAG does not violate hierarchical constraints imposed by the GO MF hierarchy in cross-validated datasets, just replace replace ``-e ho`` with ``-e cv`` in the :ref:`hemdagcheck` script:

.. code-block:: bash

    Rscript hemdag-checker.R -e cv

Evaluation
-------------------
To evaluate HEMDAG in cross-validated experiments just replace ``-e ho`` with ``-e cv`` in the :ref:`hemdageval` script:

.. code-block:: bash

    ## single evaluation call
    Rscript hemdag-eval.R -o 7227_drome -d mf -e cv -f ranger -a gpav

    ## generate chunk evaluation calls
    perl hemdag-chunk.pl -e cv -c 6 > hemdag-cv-eval.sh

    ## evaluate HEMDAG in chunks
    bash hemdag-cv-eval.sh > out


