#############
## TPR-DAG ##
#############

#' @name TPR-DAG-variants
#' @seealso \code{\link{GPAV-DAG}}, \code{\link{HTD-DAG}}
#' @title TPR-DAG Ensemble Variants
#' @description Function gathering the true-path-rule-based hierarchical learning ensemble algorithms and its variants. 
#' In their more general form the \code{TPR-DAG} algorithms adopt a two step learning strategy:
#' \enumerate{
#'  \item in the first step they compute a \emph{per-level bottom-up} visit from the leaves to the root to propagate positive predictions 
#'  across the hierarchy;
#'  \item in the second step they compute a \emph{per-level top-down} visit from the root to the leaves in order to assure the hierarchical 
#'  consistency of the predictions;
#' }
#' @details The \emph{vanilla} \code{TPR-DAG} adopts a per-level bottom-up traversal of the DAG to correct the flat predictions \eqn{\hat{y}_i}:
#' \deqn{
#'  \bar{y}_i := \frac{1}{1 + |\phi_i|} (\hat{y}_i + \sum_{j \in \phi_i} \bar{y}_j)
#' }
#' where \eqn{\phi_i} are the positive children of \eqn{i}.
#' Different strategies to select the positive children \eqn{\phi_i} can be applied:
#' \enumerate{
#'  \item \strong{Threshold-Free} strategy: the positive nodes are those children that can increment the score of the node \eqn{i}, that is those nodes 
#'  that achieve a score higher than that of their parents:
#'  \deqn{
#'      \phi_i := \{ j \in child(i) | \bar{y}_j > \hat{y}_i \}
#'  }
#'  \item \strong{Threshold} strategy: the positive children are selected on the basis of a threshold that can be selected in two different ways:
#'  \enumerate{
#'      \item for each node a constant threshold \eqn{\bar{t}} is a priori selected:
#'      \deqn{
#'          \phi_i := \{ j \in child(i) | \bar{y}_j > \bar{t} \}
#'      }
#'      For instance if the predictions represent probabilities it could be meaningful to a priori select \eqn{\bar{t}=0.5}.
#'      \item the threshold is selected to maximize some performance metric \eqn{\mathcal{M}} estimated on the training data, as for instance
#'      the Fmax or the AUPRC. In other words the threshold is selected to maximize some measure of accuracy of the predictions 
#'      \eqn{\mathcal{M}(j,t)} on the training data for the class \eqn{j} with respect to the threshold \eqn{t}. 
#'      The corresponding set of positives \eqn{\forall i \in V} is:
#'      \deqn{
#'          \phi_i := \{ j \in child(i) | \bar{y}_j > t_j^*,  t_j^* = \arg \max_{t} \mathcal{M}(j,t) \}
#'      }
#'      For instance \eqn{t_j^*} can be selected from a set of \eqn{t \in (0,1)} through internal cross-validation techniques.
#'  }
#' }
#' The weighted \code{TPR-DAG} version can be designed by adding a weight \eqn{w \in [0,1]} to balance between the 
#' contribution of the node \eqn{i} and that of its positive children \eqn{\phi}, through their convex combination:
#' \deqn{
#'  \bar{y}_i := w \hat{y}_i + \frac{(1 - w)}{|\phi_i|} \sum_{j \in \phi_i} \bar{y}_j
#' }
#' If \eqn{w=1} no weight is attributed to the children and the \code{TPR-DAG} reduces to the \code{HTD-DAG} algorithm, since in this
#' way only the prediction for node \eqn{i} is used in the bottom-up step of the algorithm. If \eqn{w=0} only the predictors 
#' associated to the children nodes vote to predict node \eqn{i}. In the intermediate cases we attribute more importance to the predictor for the
#' node \eqn{i} or to its children depending on the values of \eqn{w}.
#'
#' The contribution of the descendants of a given node decays exponentially with their distance from the node itself. To enhance the 
#' contribution of the most specific nodes to the overall decision of the ensemble we designed a novel variant that we named \code{DESCENS}. 
#' The novelty of \code{DESCENS} consists in strongly considering the contribution of all the descendants of each node instead of 
#' only that of its children. Therefore \code{DESCENS} predictions are more influenced by the information embedded in the leaves nodes, 
#' that are the classes containing the most informative and meaningful information from a biological and medical standpoint. 
#' For the choice of the ``positive'' descendants we use the same strategies adopted for the selection of the ``positive'' 
#' children shown above. Furthermore, we designed a variant specific only for \code{DESCENS}, that we named \code{DESCENS}-\eqn{\tau}.
#' The \code{DESCENS}-\eqn{\tau} variants balances the contribution between the ``positives'' children of a node \eqn{i} 
#' and that of its ``positives'' descendants excluding its children by adding a weight \eqn{\tau \in [0,1]}:
#' \deqn{
#' \bar{y}_i := \frac{\tau}{ 1 +|\phi_i|} ( \hat{y}_i + \sum_{j \in \phi_i} \bar{y}_j ) + \frac{1-\tau}{1+|\delta_i|} ( \hat{y}_i + \sum_{j\in \delta_i} \bar{y}_j )
#' }
#' where \eqn{\phi_i} are the ``positive'' children of \eqn{i} and \eqn{\delta_i=\Delta_i \setminus \phi_i} the descendants of \eqn{i} without its children. 
#' If \eqn{\tau=1} we consider only the contribution of the ``positive'' children of \eqn{i}; if \eqn{\tau=0} only the descendants that are not
#' children contribute to the score, while for intermediate values of \eqn{\tau} we can balance the contribution of \eqn{\phi_i} and 
#' \eqn{\delta_i} positive nodes.
#'
#' Simply by replacing the \code{HTD} (\code{\link{HTD-DAG}}) top-down step with the \code{GPAV} approach (\code{\link{GPAV-DAG}}) we can design the
#' \code{TPR-DAG} variant \code{ISO-TPR}. The most important feature of \code{ISO-TPR} is that it maintains the hierarchical constraints by
#' construction and selects the closest solution (in the least square sense) to the bottom-up predictions that obeys the true path rule.
#' Obviously, any aforementioned strategy for the selection of ``positive'' children or descendants can be applied before executing the \code{GPAV} correction.
#' @param S a named flat scores matrix with examples on rows and classes on columns.
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param root name of the class that it is on the top-level of the hierarchy (\code{def. root="00"}).
#' @param positive choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#'  \item \code{children} (\code{def.}): for each node are considered its positive children;
#'  \item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#'  \item \code{threshold.free} (\code{def.}): positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#'  \item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#'  \item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#'  \item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#'  \item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy. 
#'  NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positive=descendants};
#' }
#' @param topdown strategy to make the scores hierarchy-consistent. It can be one of the following values:
#' \itemize{
#'  \item \code{htd} (\code{def.}): \code{HTD-DAG} strategy is applied (\code{\link{HTD-DAG}});
#'  \item \code{gpav}: \code{GPAV} strategy is applied (\code{\link{GPAV-DAG}});
#' }
#' @param t threshold for the choice of positive nodes (\code{def. t=0}). Set \code{t} only for the variants that requiring 
#' a threshold for the selection of the positive nodes, otherwise set \code{t} to zero.
#' @param w weight to balance between the contribution of the node \eqn{i} and that of its positive nodes. Set \code{w} only for the
#' \emph{weighted} variants, otherwise set \code{w} to zero.
#' @param W vector of weight relative to a single example. If \code{W=NULL} (def.) it is assumed that \code{W} is a unitary vector of the same length of 
#' the columns' number of the matrix \code{S} (root node included). Set \code{W} only if \code{topdown=GPAV}.
#' @param parallel boolean value:
#' \itemize{
#'  \item \code{TRUE}: execute the parallel implementation of GPAV (\code{\link{GPAV.parallel}});
#'  \item \code{FALSE} (def.): execute the sequential implementation of GPAV (\code{\link{GPAV.over.examples}});
#' }
#' Use \code{parallel} only if \code{topdown=GPAV}; otherwise set \code{parallel=FALSE}.
#' @param ncores number of cores to use for parallel execution (\code{def. 8}). Set \code{ncores=1} if \code{parallel=FALSE}, 
#' otherwise set \code{ncores} to the desired number of cores.
#' Use \code{ncores} if and only if \code{topdown=GPAV}; otherwise set \code{parallel=1}.
#' @return a named matrix with the scores of the classes corrected according to the TPR-DAG ensemble algorithm.
#' @export 
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' root <- root.node(g);
#' S.tpr <- tpr.dag(S, g, root, positive="children", bottomup="threshold.free", topdown="htd", t=0, w=0, W=NULL, parallel=FALSE, ncores=1);
tpr.dag <- function(S, g, root="00", positive="children", bottomup="threshold.free", topdown="htd", t=0, w=0, W=NULL, parallel=FALSE, ncores=1){
    ## parameters check
    if(positive!="children" && positive!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && 
        bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau" || topdown!="HTD" && topdown!="GPAV")
        stop("tpr.dag: positive or bottomup or topdown value misspelled", call.=FALSE);
    if(positive=="children" && bottomup=="tau")
        stop("tpr.dag: tau is a descendants variants. Please set positive to descendants", call.=FALSE);
    if(bottomup=="threshold" || bottomup=="tau")
        w <- 0;
    if(bottomup=="threshold.free"){
        t <- 0;
        w <- 0;
    }
    if(bottomup=="weighted.threshold.free")
        t <-0;
    if(t==1 || w==1)
        warning("tpr.dag: when t or w is equal to 1, TPR-DAG is reduced to HTD-DAG", call.=FALSE);  
    if(topdown=="gpav" && parallel==TRUE && ncores<2)
        warning("gpav: set ncores greater than 2 to exploit the gpav parallel version", call.=FALSE);
    if(topdown=="gpav" && parallel==FALSE && ncores>=2)
        warning("tpr.dag: set ncores greater than 2 to exploit the gpav parallel version", call.=FALSE);
    if(topdown=="htd" && (parallel==TRUE || ncores>=2))
        warning("TPR-DAG: does not exist a parallel version of HTD. Set 'parallel' to FALSE and/or 'ncores' to 1 to avoid this warning message", call.=FALSE);   

    ## add root node to S if it does not exist
    if(!(root %in% colnames(S))){
        max.score <- max(S);
        z <- rep(max.score,nrow(S));
        S <- cbind(z,S);
        colnames(S)[1] <- root;
    }
    ## check consistency between nodes of g and classes of S
    class.check <- ncol(S)!=numNodes(g);
    if(class.check)
        stop("tpr.dag: mismatch between the number of nodes of the graph and the number of class of the scores matrix", call.=FALSE);

    ## computing graph levels
    levels <- graph.levels(g,root);

    ## bottom-up visit: positive children selection
    if(positive=="children"){
        chd.bup <- get.children.bottom.up(g,levels);
        for(i in 1:length(chd.bup)){
            if(length(chd.bup[[i]])!=0){
                parent <- S[,names(chd.bup[i])];
                children <- as.matrix(S[,chd.bup[[i]]]);
                # colnames(children) <- chd.bup[[i]]
                for(j in 1:length(parent)){
                    if(bottomup=="threshold"){
                        child.set <- children[j,] > t;    # positive children selection
                        child.pos <- children[j,][child.set];
                        parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat scores correction
                    }else if(bottomup=="threshold.free"){
                        child.set <- children[j,] > parent[j]; # positive children selection
                        child.pos <- children[j,][child.set];
                        parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat score correction
                    }else if(bottomup=="weighted.threshold.free"){
                        child.set <- children[j,] > parent[j];    # positive children selection
                        child.pos <- children[j,][child.set];
                        if(length(child.pos)!=0){
                            parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score correction
                        }
                    }else if(bottomup=="weighted.threshold"){
                        child.set <- children[j,] > t;    # positive children selection
                        child.pos <- children[j,][child.set];
                        if(length(child.pos)!=0){
                            parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score prediction
                        }
                    }
                }   
                S[,names(chd.bup[i])] <- parent;
            }
        }
    ## bottom-up visit: positive descendants selection
    }else if(positive=="descendants"){
        if(bottomup=="tau")
            chd.bup <- get.children.bottom.up(g,levels);
        desc.bup <- build.descendants.bottom.up(g,levels);
        nodes <- names(desc.bup);
        for(i in 1:length(desc.bup)){
            if(length(desc.bup[[i]])!=1){
                node.curr <- nodes[i];
                parent <- S[,names(desc.bup[i])];
                tmp <- setdiff(desc.bup[[i]],node.curr);
                if(bottomup=="tau"){
                    delta <- setdiff(tmp, chd.bup[[i]]);  # descendants without children 
                    children <- as.matrix(S[,chd.bup[[i]]]);    # genes considering children node 
                    desc <-  as.matrix(S[,delta]);      # genes considering descendants nodes without children
                }else{
                    desc <- as.matrix(S[,tmp]);
                }
                for(j in 1:length(parent)){
                    if(bottomup=="threshold"){
                        desc.set <- desc[j,] > t;    # positive descendants selection
                        desc.pos <- desc[j,][desc.set];
                        parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
                    }else if(bottomup=="threshold.free"){
                        desc.set <- desc[j,] > parent[j];   # positive descendants selection
                        desc.pos <- desc[j,][desc.set];
                        parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
                    }else if(bottomup=="weighted.threshold.free"){
                        desc.set <- desc[j,] > parent[j];
                        desc.pos <- desc[j,][desc.set];
                        if(length(desc.pos)!=0){
                            parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
                        }
                    }else if(bottomup=="weighted.threshold"){
                        desc.set <- desc[j,] > t;
                        desc.pos <- desc[j,][desc.set];
                        if(length(desc.pos)!=0){
                            parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
                        }
                    }else if(bottomup=="tau"){
                        desc.set <- desc[j,] > parent[j];           # positive descendants (without children) selection
                        desc.pos <- desc[j,][desc.set];
                        child.set <- children[j,] > parent[j];      # positive children selection
                        child.pos <- children[j,][child.set];
                        parent[j] <- t * ((parent[j] + sum(child.pos))/(1+length(child.pos))) + (1-t) * ((parent[j] + sum(desc.pos))/(1+length(desc.pos)));
                    }
                    S[,names(desc.bup[i])] <- parent;
                }
            }
        }
    }
    # top-down visit
    if(topdown=="htd"){
        par.tod <- get.parents.top.down(g,levels,root);
        for(i in 1:length(par.tod)){
            child <- S[,names(par.tod[i])];
            parents <- as.matrix(S[,par.tod[[i]]]);
            # colnames(parents) <- par.tod[[i]]
            # Note: the version with an apply and an ifelse statement is slower ...
            for(j in 1:length(child)){
                x <- min(parents[j,]);
                if(x < child[j]){
                    child[j] <- x;    # hierarchical correction
                }
            }
            S[,names(par.tod[i])] <- child;
        }
    }else if(topdown=="gpav"){
        if(parallel){
            S <- gpav.parallel(S, g, W=W, ncores=ncores);
        }else{
            S <- gpav.over.examples(S, W=W, g);
        }
    } 
    return(S);
}

#' @name TPR-DAG-cross-validation
#' @title TPR-DAG cross-validation experiments
#' @seealso \code{\link{TPR-DAG-variants}}
#' @description Correct the computed scores in a hierarchy according to the chosen TPR-DAG ensemble variant. 
#' @details The parametric hierarchical ensemble variants are cross-validated maximizing the parameter on the metric selected in \code{metric}, 
#' @param S a named flat scores matrix with examples on rows and classes on columns.
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param norm boolean value. Should the flat score matrix be normalized? By default \code{norm=FALSE}. If \code{norm=TRUE} the matrix \code{S} is normalized according to \code{norm.type}.
#' @param norm.type can be one of the following values: 
#'  \enumerate{
#'  \item \code{NULL} (def.): none normalization is applied (\code{norm=FALSE})
#'  \item \code{maxnorm}: each score is divided for the maximum value of each class;
#'  \item \code{qnorm}: quantile normalization. \pkg{preprocessCore} package is used; 
#'  }
#' @param W vector of weight relative to a single example. If \code{W=NULL} (def.) it is assumed that \code{W} is a unitary vector of the same length of 
#' the columns' number of the matrix \code{S} (root node included). Set \code{W} only if \code{topdown=GPAV}.
#' @param parallel boolean value. Should the parallel version \code{GPAV-DAG} be run?
#' \itemize{
#'    \item \code{TRUE}: execute the parallel implementation of \code{GPAV} (\code{\link{gpav.parallel}});
#'    \item \code{FALSE} (\code{def.}): execute the sequential implementation of \code{GPAV} (\code{\link{gpav.over.examples}});
#' }
#' @param ncores number of cores to use for parallel execution (\code{def. 1}).
#' @param positive choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#'  \item \code{children} (\code{def.}): for each node are considered its positive children;
#'  \item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#'  \item \code{threshold.free} (\code{def.}): positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#'  \item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#'  \item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#'  \item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#'  \item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy;
#'  NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positive=descendants};
#' }
#' @param topdown strategy to make the scores 'ontology-aware'. It can be one of the following values:
#' \itemize{
#'  \item \code{htd} (\code{def.}): \code{HTD-DAG} strategy is applied (\code{\link{HTD-DAG}});
#'  \item \code{gpav}: \code{GPAV} strategy is applied (\code{\link{GPAV-DAG}});
#' }
#' @param threshold range of threshold values to be tested in order to find the best threshold (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but the execution time will be higher.
#' Set this parameter only for the \emph{thresholded} variants; for the \emph{threshold-free} variants, \code{threshold} is automatically set to zero.
#' @param weight range of weight values to be tested in order to find the best weight (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but obviously the execution time will be higher.
#' Set this parameter only for the \emph{weighted} variants; for the \emph{weight-free} variants,
#' the parameter \code{weight} is automatically set to zero.
#' @param metric a string character specifying the performance metric on which maximizing the parametric ensemble variant. 
#' It can be one of the following values:
#' \enumerate{
#' \item \code{prc}: the parametric ensemble variant is maximized on the basis of AUPRC (\code{\link{AUPRC}});
#' \item \code{fmax}: the parametric ensemble variant is maximized on the basis of Fmax (\code{\link{Multilabel.F.measure}};
#' \item \code{NULL}: on the \code{threshold.free} variant none parameter optimization is needed, since the variant is non-parametric.
#' So, if \code{bottomup=threshold.free} set \code{metric=NULL} (\code{def.});
#' }
#' @param kk number of folds of the cross validation (\code{def: kk=5}) on which tuning the parameters \code{threshold}, \code{weight} and 
#' \code{tau} of the parametric ensemble variants. For the non-parametric variants(i.e. if \code{bottomup = threshold.free}), \code{kk} is automatically set to zero. 
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. If \code{compute.performance=FALSE} and \code{bottomup=threshold.free}, \code{seed} is automatically set to \code{NULL}.
#' @return a named matrix with the scores of the classes corrected according to the chosen TPR-DAG ensemble algorithm.
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' S.tpr <- tpr.dag.cv(S, g, norm=FALSE, norm.type=NULL, parallel=FALSE, ncores=1, W=NULL, positive="children", bottomup="threshold.free", topdown="htd", threshold=0, weight=0, W=NULL, metric=NULL);
tpr.dag.cv <- function(S, g, norm=FALSE, norm.type=NULL, parallel=FALSE, ncores=1, W=NULL, positive="children", bottomup="threshold.free", topdown="HTD",  
    threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=23, metric=NULL){
    ## parameters check 
    if(positive!="children" && positive!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau" || topdown!="htd" && topdown!="gpav")
        stop("tpr.dag.cv: positive or bottomup or topdown value misspelled", call.=FALSE);
    if(positive=="children" && bottomup=="tau")
        stop("tpr.dag.cv: tau is a descendants variants. Please set positive to descendants", call.=FALSE);
    if(bottomup=="threshold" || bottomup=="tau")
        weight <- 0;
    if(bottomup=="threshold.free"){
        threshold <- 0; 
        weight <- 0;
    }
    if(bottomup=="weighted.threshold.free")
        threshold <- 0;
    if(norm==TRUE && is.null(norm.type))
        stop("tpr.dag.cv: choose a normalization methods among those available", call.=FALSE);
    if(norm==FALSE && !is.null(norm.type))
        warning("tpr.dag.cv: ", paste0("set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);   
    if((is.null(kk) || kk<=1) && bottomup!="threshold.free")
        stop("tpr.dag.cv: smallest number of folds to define test and training set is 2. Set kk larger or equal to 2", call.=FALSE);
    if(!is.null(kk) && bottomup=="threshold.free")
        kk <- NULL; 
    if(metric!="fmax" && metric!="prc" && !is.null(metric))
        stop("tpr.dag.cv: value of parameter 'metric' misspelled", call.=FALSE);
    if(is.null(metric) && bottomup!="threshold.free")
        stop(paste0("tpr.dag.cv: the bottom-up approach '", bottomup, "' is parametric"),". Select the metric on which maximize according to those available", call.=FALSE); 
    if(!is.null(metric) && bottomup=="threshold.free") 
        warning("tpr.dag.cv: the bottom-up approach 'threshold.free' is non-parametric. Set metric to NULL to avoid this warning message", call.=FALSE);
    if(is.null(seed) && bottomup!="threshold.free")
        warning("tpr.dag.cv: folds are generate without seed initialization", call.=FALSE);
    
    ## normalization
    if(norm){
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "normalization: done", "\n");
    }
    
    ## compute root node
    root <- root.node(g);

    ## tpr-dag hierarchical correction 
    if(bottomup=="threshold.free"){
        S.hier <- tpr.dag(S, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, t=0, w=0, W=W, parallel=parallel, ncores=ncores);
        cat("tpr-dag correction: done\n");
        rm(S); gc();
    }else{
        ## let's start k-fold crossing validation for choosing best threshold and weight maximizing on the selected metric
        testIndex <- do.unstratified.cv.data(S, kk=kk, seed=seed); 
        S.hier <- c(); # variable to host the k-assembled sub-matrix  
        # training.top <- vector(mode="list", length=kk); ## for check
        for(k in 1:kk){
            ## training test
            training <- S[-testIndex[[k]],];
            target.training <- ann[-testIndex[[k]],];           
            ## test set
            test <- S[testIndex[[k]],]; 
            target.test <- ann[testIndex[[k]],];
            ## metric initialization        
            top.metric <- 0;
            bestT <- 0;
            bestW <- 0;
            for(t in threshold){
                for(w in weight){
                    pred.training <- TPR.DAG(training, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, w=w, t=t, W=W, parallel=parallel, ncores=ncores);
                    if(metric=="fmax"){
                        if(root %in% colnames(pred.training))
                            pred.training <- pred.training[,-which(colnames(pred.training)==root)];
                        training.metric <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE)[["F"]];
                    }else{
                        if(root %in% colnames(pred.training))
                            pred.training <- pred.training[,-which(colnames(pred.training)==root)];
                        training.metric <- auprc.single.over.class(target.training, pred.training, folds=NULL, seed=NULL)$average;
                    }
                    if(training.metric > top.metric){
                        top.metric <- training.metric;
                        bestT <- t;
                        bestW <- w;
                        # training.top[[k]] <- c(metric=top.metric, best.thres=bestT, best.weight=bestW);
                        if(bottomup=="threshold" || bottomup=="tau"){
                            cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best threshold:", bestT, sep="\t", "\n");
                        }else if(bottomup=="weighted.threshold.free"){
                            cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best weight:", bestW, sep="\t", "\n");
                        }
                        else{
                            cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best threshold:",bestT, "best weight:", bestW, sep="\t", "\n");
                        }
                    }
                }
            }
            ## test set 
            pred.test <- tpr.dag(test, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, t=bestT, w=bestW, W=W, parallel=parallel, ncores=ncores);
            ## assembling the hierarchical scores of each k sub-matrix
            S.hier <- rbind(S.hier, pred.test); 
        }
        ## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the beginning matrix
        S.hier <- S.hier[rownames(S),];
        cat("tpr-dag correction: done\n");
        rm(S, testIndex, pred.test, test, training, target.test, target.training); gc();
    }
    return(S.hier);
}

#' @name TPR-DAG-holdout
#' @title TPR-DAG holdout experiments
#' @seealso \code{\link{TPR-DAG-variants}}
#' @description Correct the computed scores in a hierarchy according to the chosen TPR-DAG ensemble variant by applying a classical holdout procedure.
#' @details The parametric hierarchical ensemble variants are cross-validated maximizing the parameter on the metric selected in \code{metric}, 
#' @param S a named flat scores matrix with examples on rows and classes on columns.
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param testIndex a vector of integer numbers corresponding to the indexes of the elements (rows) of the scores matrix \code{S} to be used in the test set.
#' @param norm boolean value. Should the flat score matrix be normalized? By default \code{norm=FALSE}. If \code{norm=TRUE} the matrix \code{S} is normalized according to \code{norm.type}.
#' @param norm.type can be one of the following values: 
#'  \enumerate{
#'  \item \code{NULL} (def.): none normalization is applied (\code{norm=FALSE})
#'  \item \code{maxnorm}: each score is divided for the maximum value of each class;
#'  \item \code{qnorm}: quantile normalization. \pkg{preprocessCore} package is used; 
#'  }
#' @param W vector of weight relative to a single example. If \code{W=NULL} (def.) it is assumed that
#' \code{W} is a unitary vector of the same length of the columns' number of the matrix \code{S} (root node included).
#' @param parallel boolean value. Should the parallel version \code{GPAV-DAG} be run?
#' \itemize{
#'    \item \code{TRUE}: execute the parallel implementation of \code{GPAV} (\code{\link{gpav.parallel}});
#'    \item \code{FALSE} (\code{def.}): execute the sequential implementation of \code{GPAV} (\code{\link{gpav.over.examples}});
#' }
#' @param ncores number of cores to use for parallel execution (\code{def. 8}).
#' @param positive choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#'  \item \code{children} (\code{def.}): for each node are considered its positive children;
#'  \item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#'  \item \code{threshold.free} (\code{def.}): positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#'  \item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#'  \item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#'  \item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#'  \item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy;
#'  NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positive=descendants};
#' }
#' @param topdown strategy to make the scores 'ontology-aware'. It can be one of the following values:
#' \itemize{
#'  \item \code{htd} (\code{def.}): \code{HTD-DAG} strategy is applied (\code{\link{HTD-DAG}});
#'  \item \code{gpav}: \code{GPAV} strategy is applied (\code{\link{GPAV-DAG}});
#' }
#' @param threshold range of threshold values to be tested in order to find the best threshold (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but the execution time will be higher.
#' Set this parameter only for the \emph{thresholded} variants; for the \emph{threshold-free} variants, \code{threshold} is automatically set to zero.
#' @param weight range of weight values to be tested in order to find the best weight (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but obviously the execution time will be higher.
#' Set this parameter only for the \emph{weighted} variants; for the \emph{weight-free} variants,
#' the parameter \code{weight} is automatically set to zero.
#' @param metric a string character specifying the performance metric on which maximizing the parametric ensemble variant. 
#' It can be one of the following values:
#' \enumerate{
#' \item \code{prc}: the parametric ensemble variant is maximized on the basis of AUPRC (\code{\link{AUPRC}});
#' \item \code{fmax}: the parametric ensemble variant is maximized on the basis of Fmax (\code{\link{Multilabel.F.measure}};
#' \item \code{NULL}: on the \code{threshold.free} variant none parameter optimization is needed, since the variant is non-parametric.
#' So, if \code{bottomup=threshold.free} set \code{metric=NULL} (\code{def.});
#' }
#' @param kk number of folds of the cross validation (\code{def: kk=5}) on which tuning the parameters \code{threshold}, \code{weight} and 
#' \code{tau} of the parametric ensemble variants. For the non-parametric variants(i.e. if \code{bottomup = threshold.free}), \code{kk} is automatically set to zero. 
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. If \code{compute.performance=FALSE} and \code{bottomup=threshold.free}, \code{seed} is automatically set to \code{NULL}.
#' @return a named matrix with the scores of the classes corrected according to the chosen TPR-DAG ensemble algorithm. Rows of the matrix are shrunk to \code{testIndex}.
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(test.index);
#' S.tpr <- tpr.dag.holdout(S, g, testIndex=test.index, norm=FALSE, norm.type=NULL, parallel=FALSE, ncores=1, W=NULL, 
#' positive="children", bottomup="threshold.free", topdown="htd", threshold=0, weight=0, W=NULL, metric=NULL);
tpr.dag.holdout <- function(S, g, testIndex, norm=FALSE, norm.type=NULL, W=NULL, parallel=FALSE, ncores=1, positive="children", bottomup="threshold.free",
    topdown="htd",  threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=23, metric=NULL){
    ## parameters check 
    if(positive!="children" && positive!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau" || topdown!="htd" && topdown!="gpav")
        stop("tpr.dag.cv: positive or bottomup or topdown value misspelled", call.=FALSE);
    if(positive=="children" && bottomup=="tau")
        stop("tpr.dag.cv: tau is a descendants variants. Please set positive to descendants", call.=FALSE);
    if(bottomup=="threshold" || bottomup=="tau")
        weight <- 0;
    if(bottomup=="threshold.free"){
        threshold <- 0; 
        weight <- 0;
    }
    if(bottomup=="weighted.threshold.free")
        threshold <- 0;
    if(norm==TRUE && is.null(norm.type))
        stop("tpr.dag.cv: choose a normalization methods among those available", call.=FALSE);
    if(norm==FALSE && !is.null(norm.type))
        warning("tpr.dag.cv: ", paste0("set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);   
    if((is.null(kk) || kk<=1) && bottomup!="threshold.free")
        stop("tpr.dag.cv: smallest number of folds to define test and training set is 2. Set kk larger or equal to 2", call.=FALSE);
    if(!is.null(kk) && bottomup=="threshold.free")
        kk <- NULL; 
    if(metric!="fmax" && metric!="prc" && !is.null(metric))
        stop("tpr.dag.cv: value of parameter 'metric' misspelled", call.=FALSE);
    if(is.null(metric) && bottomup!="threshold.free")
        stop(paste0("tpr.dag.cv: the bottom-up approach '", bottomup, "' is parametric"),". Select the metric on which maximize according to those available", call.=FALSE); 
    if(!is.null(metric) && bottomup=="threshold.free") 
        warning("tpr.dag.cv: the bottom-up approach 'threshold.free' is non-parametric. Set metric to NULL to avoid this warning message", call.=FALSE);
    if(is.null(seed) && bottomup!="threshold.free")
        warning("tpr.dag.cv: folds are generate without seed initialization", call.=FALSE);

    ## normalization
    if(norm){
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "normalization: done", "\n");
    }
    
    ## compute root node
    root <- root.node(g);

    ## check if |flat matrix classes| = |annotation matrix classes| 
    ## if not the classes of annotation matrix are shrunk to those of flat matrix
    # class.check <- ncol(S)!=ncol(ann);
    # if(class.check){
    #     ann <- ann[,colnames(S)];
    #     nd <- c(root, colnames(S));
    #     g <- do.subgraph(nd, g, edgemode="directed");
    # }

    ## scores flat matrix are shrunk to test and training test respectively
    S.test <- S[testIndex,];
    S.training <- S[-testIndex,];

    ## tpr-dag hierarchical correction 
    if(bottomup=="threshold.free"){
        S.hier <- tpr.dag(S.test, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, t=0, w=0, W=W, parallel=parallel, ncores=ncores);
        cat("tpr-dag correction: done\n");
        rm(S);
    }else{
        ## let's start k-fold crossing validation for choosing best threshold and weight maximizing on the selected metric
        foldIndex <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
        # training.top <- vector(mode="list", length=kk); ## for check
        for(k in 1:kk){
            ## training and test set
            training <- S.training[foldIndex[[k]],];
            target.training <- ann.training[foldIndex[[k]],];
            ## metric initialization        
            top.metric <- 0;
            bestT <- 0;
            bestW <- 0;
            for(t in threshold){
                for(w in weight){                    
                    pred.training <- tpr.dag(training, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, w=w, t=t, W=W, parallel=parallel, ncores=ncores);
                    if(metric=="FMAX"){
                        if(root %in% colnames(pred.training))
                            pred.training <- pred.training[,-which(colnames(pred.training)==root)];
                        training.metric <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE)[["F"]];
                    }else{
                        if(root %in% colnames(pred.training))
                            pred.training <- pred.training[,-which(colnames(pred.training)==root)];
                        training.metric <- auprc.single.over.class(target.training, pred.training, folds=NULL, seed=NULL)$average;     
                    }
                    if(training.metric > top.metric){
                        top.metric <- training.metric;
                        bestT <- t;
                        bestW <- w;
                        # training.top[[k]] <- c(metric=top.metric, best.thres=bestT, best.weight=bestW);
                        if(bottomup=="threshold" || bottomup=="tau"){
                            cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best threshold:", bestT, sep="\t", "\n");
                        }else if(bottomup=="weighted.threshold.free"){
                            cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best weight:", bestW, sep="\t", "\n");
                        }
                        else{
                            cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best threshold:",bestT, "best weight:", bestW, sep="\t", "\n");
                        }
                    }
                }
            }
        }
        S.hier <- tpr.dag(S.test, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, t=bestT, w=bestW, W=W, parallel=parallel, ncores=ncores);
        cat("tpr-dag correction: done\n");
        rm(S, S.test, S.training, training); gc();
    }
    return(S.hier);
}

#' @title Unstratified Cross Validation
#' @description This function splits a dataset in k-fold in an unstratified way, i.e. a fold does not contain an equal amount of positive and 
#' negative examples. This function is used to perform k-fold cross-validation experiments in a hierarchical correction contest where 
#' splitting dataset in a stratified way is not needed. 
#' @param S matrix of the flat scores. It must be a named matrix, where rows are example (e.g. genes) and columns are classes/terms (e.g. GO terms).
#' @param kk number of folds in which to split the dataset (\code{def. k=5}).
#' @param seed seed for the random generator. If \code{NULL} (def.) no initialization is performed.
#' @return a list with \eqn{k=kk} components (folds). Each component of the list is a character vector contains the index of the examples, i.e. the 
#' index of the rows of the matrix S.
#' @export
#' @examples
#' data(scores);
#' foldIndex <- do.unstratified.cv.data(S, kk=5, seed=23);
do.unstratified.cv.data <- function(S, kk=5, seed=NULL){
    set.seed(seed);
    examples <- 1:nrow(S);
    n <- nrow(S);
    size <- c();
    folds <- vector(mode="list", length=kk)
    names(folds) <- paste0(rep("fold",kk), 1:kk)
    for(k in 1:kk){
        first <- ((k - 1) * n) %/% kk
        last <- (k * n) %/% kk
        size <- last-first;
        x    <- sample(examples,size);
        folds[[k]] <- x;
        examples <- setdiff(examples,x);
    }
    return(folds);
}
