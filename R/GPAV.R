##############
## GPAV-DAG ##
##############

#' @title Binary Upper Triangular Adjacency Matrix
#' @description This function returns a binary square upper triangular matrix where rows and columns correspond to the nodes' name of the graph \code{g}.
#' @details The nodes of the matrix are topologically sorted (by using the \code{tsort} function of the \pkg{RBGL} package). 
#' Let's denote with \code{adj} our adjacency matrix. Then \code{adj} represents a partial
#' order data set in which the class \code{j} dominates the class \code{i}. In other words, \code{adj[i,j]=1} means that \code{j} dominates \code{i};
#' \code{adj[i,j]=0} means that there is no edge between the class \code{i} and the class \code{j}. Moreover the nodes of \code{adj} are
#' ordered such that \code{adj[i,j]=1} implies \eqn{i < j}, i.e. \code{adj} is upper triangular.
#' @param g a graph of class \code{graphNELL} representing the hierarchy of the class.
#' @return an adjacency matrix which is square, logical and upper triangular.
#' @export
#' @examples
#' data(graph);
#' adj <- adj.upper.tri(g);
adj.upper.tri <- function(g){
    ## 0. write the graph g as pair parent-child (source- destination)
    num.ed <- numEdges(g);
    num.nd <- numNodes(g);
    m <- matrix(character(num.ed*2), ncol=2);
    ed <- edges(g);
    count <- 0;
    par <- names(ed);
    for(i in 1:num.nd){
        children <- ed[[i]];
        len.x <- length(children);
        if(len.x!=0){
            for(j in 1:len.x){
                count <- count + 1;
                m[count,] <- c(par[i], children[j]);
            }
        }
    }

    ## 1. topological sorting: nodes are ordering in a linear way such as vertex u comes before v.
    tsort.nd <- tsort(g);

    ## 2. map each node of the graph to a integer number (in a decreasing order)
    source <- mapvalues(m[,1], from=tsort.nd[1:length(tsort.nd)], to=length(tsort.nd):1, warn_missing=F);
    destination  <- mapvalues(m[,2], from=tsort.nd[1:length(tsort.nd)], to=length(tsort.nd):1, warn_missing=F);

    ## 3. build upper triangular logical adjacency constraints matrix. this matrix should be sparse
    eM <- cbind(from=as.numeric(destination), to=as.numeric(source));
    adj <- matrix(0, nrow=num.nd, ncol=num.nd);
    rev.nd <- tsort.nd[num.nd:1];  ## reverse topological sorting for a right mapping with number
    dimnames(adj) <- list(rev.nd, rev.nd);
    adj[eM] <- 1;
    return(adj);
}

#' @title Generalized Pool-Adjacent Violators
#' @description Implementation of \code{GPAV} (Generalized Pool-Adjacent Violators) algorithm.
#' (\cite{Burdakov et al., In: Di Pillo G, Roma M, editors. An O(n2) Algorithm for Isotonic Regression. Boston, MA: Springer US; 2006. 
#' p. 25–33. Available from: \href{https://doi.org/10.1007/0-387-30065-1_3}{https://doi.org/10.1007/0-387-30065-1_3}} 
#' @details Given the constraints adjacency matrix of the graph, a vector of scores \eqn{\hat{y} \in R^n} and a vector of strictly positive
#' weights \eqn{w \in R^n}, the \code{GPAV} algorithm returns a vector \eqn{\bar{y}} which is as close as possible, in the least-squares sense,
#' to the response vector \eqn{\hat{y}} and whose components are partially ordered in accordance with the constraints matrix \code{adj}.
#' In other words, \code{GPAV} solves the following problem:
#' \deqn{
#'    \bar{y} = \left\{
#'    \begin{array}{l}
#'        \min \sum_{i \in N} (\hat{y}_i - \bar{y}_i )^2\\\\
#'        \forall i, \quad  j \in par(i) \Rightarrow  \bar{y}_j  \geq \bar{y}_i
#'    \end{array}
#' \right.
#'}
#' @param Y vector of scores relative to a single example. \code{Y} must be a numeric named vector, where names
#' correspond to classes' names, i.e. nodes of the graph \code{g} (root node included).
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (\code{def. W=NULL}), it is assumed that
#' \code{W} is a unitary vector of the same length of the vector \code{Y}.
#' @param adj adjacency matrix of the graph which must be sparse, logical and upper triangular. Number of columns of \code{adj} must be
#' equal to the length of \code{Y} and \code{W}.
#' @seealso \code{\link{adj.upper.tri}}
#' @return a list of 3 elements:
#' \itemize{
#'    \item \code{YFit}: a named vector with the scores of the classes corrected according to the \code{GPAV} algorithm.
#'    \code{NOTE}: the classes of \code{YFit} are topologically sorted, that is are in the same order of those of \code{adj}.
#'    \item \code{blocks}: list of vectors, containing the partitioning of nodes (represented with an integer number) into blocks;
#'    \item \code{W}: vector of weights.
#' }
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' Y <- S[3,];
#' adj <- adj.upper.tri(g);
#' S.GPAV <- GPAV(Y,W=NULL,adj);
GPAV <- function(Y, W=NULL, adj){
    if(is.null(W))
        W <- rep(1,ncol(adj));
    nW <- length(W);
    nY <- length(Y);
    nadj <- ncol(adj);
    if(nY!=nadj)
        stop("GPAV: mismatch between the number of classes between 'Y' and 'adj'", call.=FALSE);
    if(nW!=nadj)
        stop("GPAV: mismatch between the number of classes between 'Y' and 'adj'", call.=FALSE);
    ## sort nodes in the same order of adj matrix (i.e. in a topologically ordered)
    Y <- Y[colnames(adj)];
    ## just for clearnnes: we play with index and not with name...
    Y <- unname(Y);
    W <- unname(W);
    ## assign to each term an integer number, i.e. index's term
    N <- ncol(adj);
    corr <- 1:N;
    gpavres <- .C("GPAV_cpp", as.double(W), as.double(adj), as.integer(N), as.double(Y), as.integer(corr));
    YFit <- gpavres[[4]];
    corr <- gpavres[[5]];
    YFit <- YFit[corr];
    names(YFit) <- colnames(adj);
    return(YFit);
}

#' @title GPAV Over Examples
#' @description Function to compute \code{GPAV} across all the examples.
#' @param S a named flat scores matrix with examples on rows and classes on columns (root node included).
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (\code{def. W=NULL}), it is assumed that
#' \code{W} is a unitary vector of the same length of the columns' number of the matrix \code{S} (root node included).
#' @return a named matrix with the scores of the classes corrected according to the \code{GPAV} algorithm.
#' @seealso \code{\link{GPAV.parallel}}
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' S.GPAV <- GPAV.over.examples(S,W=NULL,g);
GPAV.over.examples <- function(S, g, W=NULL){
    ## check consistency between nodes of g and classes of S
    class.check <- ncol(S)!=numNodes(g);
    if(class.check)
        stop("GPAV: mismatch between the number of nodes of the graph and the number of classes of the scores matrix", call.=FALSE);
    adj <- adj.upper.tri(g);
    M <- c();
    for(i in 1:nrow(S))
        M <- rbind(M, GPAV(S[i,], W=W, adj));
    rownames(M) <- rownames(S);
    M <- M[,colnames(S)];
    S <- M; rm(M);
    return(S);
}

#' @title GPAV Over Examples -- Parallel Implementation
#' @description Function to compute \code{GPAV} across all the examples (parallel implementation).
#' @param S a named flat scores matrix with examples on rows and classes on columns (root node included).
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (\code{def. W=NULL}), it is assumed that
#' \code{W} is is a unitary vector of the same length of the columns' number of the matrix \code{S} (root node included).
#' @param ncores number of cores to use for parallel execution (\code{def. 8}).
#' @return a named matrix with the scores of the classes corrected according to the \code{GPAV} algorithm.
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' if (Sys.info()['sysname']!="Windows"){
#'    S.GPAV <- GPAV.parallel(S,W=NULL,g,ncores=2);
#' }
GPAV.parallel <- function(S, g, W=NULL, ncores=8){
    ## check consistency between nodes of g and classes of S
    class.check <- ncol(S)!=numNodes(g);
    if(class.check)
        stop("GPAV: mismatch between the number of nodes of the graph and the number of classes of the scores matrix", call.=FALSE);
    prnames <- rownames(S);
    clnames <- colnames(S);
    adj <- adj.upper.tri(g);
    if(ncores == 0){
        n.cores <- detectCores();
        if(n.cores > 3)
            ncores <- n.cores - 1;
    }
    registerDoParallel(cores=ncores);
    res.list <- foreach(i=1:nrow(S), .inorder=FALSE) %dopar% {
        res <- GPAV(S[i,], W=W, adj);
        list(protein=prnames[i], scores=res);
    }
    prnames <- unlist(lapply(res.list, '[[', 1));
    hierscores <- lapply(res.list, '[[', 2);
    S <- do.call(rbind, hierscores);
    rownames(S) <- prnames;
    S <- S[,clnames];
    rm(res.list);
    return(S);
}

#' @title GPAV -- High Level Function
#' @description High level function to correct the computed scores in a hierarchy according to the \code{GPAV} algorithm.
#' @details The function checks if the number of classes between the flat scores matrix and the annotations matrix mismatched.
#' If so, the number of terms of the annotations matrix is shrunk to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well. N.B.: it is supposed that all the nodes of the subgraph are accessible from the root.
#' @param norm boolean value:
#' \itemize{
#' \item \code{TRUE} (\code{def.}): the flat scores matrix has been already normalized in according to a normalization method;
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied;
#' }
#' @param norm.type can be one of the following three values:
#'  \enumerate{
#'  \item \code{NULL} (\code{def.}): set \code{norm.type} to \code{NULL} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used;
#'  }
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (\code{def. W=NULL}), \code{W} is a unitary
#' vector of the same length of the columns' number of the flat scores matrix (root node included).
#' @param parallel boolean value:
#' \itemize{
#'    \item \code{TRUE}: execute the parallel implementation of \code{GPAV} (\code{\link{GPAV.parallel}});
#'    \item \code{FALSE} (\code{def.}): execute the sequential implementation of \code{GPAV} (\code{\link{GPAV.over.examples}});
#' }
#' @param ncores number of cores to use for parallel execution (\code{def. 8}). Set the parameter \code{ncores} to \code{1} if the
#' parameter \code{parallel} is set to \code{FALSE}, otherwise set the desired number of cores.
#' @param folds number of folds of the cross validation on which computing the performance metrics averaged across folds (\code{def. 5}).
#' If \code{folds=NULL}, the performance metrics are computed one-shot, otherwise the performance metrics are averaged across folds.
#' If \code{compute.performance} is set to \code{FALSE}, \code{folds} is automatically set to \code{NULL}.
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. The parameter \code{seed} controls both the parameter \code{kk} and the parameter \code{folds}.
#' If \code{compute.performance} is set to \code{FALSE} and \code{bottomup} is set to \code{threshold.free}, then 
#' \code{seed} is automatically set to \code{NULL}.
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure.
#' If \code{compute.performance} is set to \code{FALSE} and \code{bottomup} is set to \code{threshold.free}, then 
#' \code{n.round} is automatically set to \code{NULL}.
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (def.): corresponds to the harmonic mean between the average precision and recall;
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples;
#' }
#' If \code{compute.performance} is set to \code{FALSE} and \code{bottomup} is set to \code{threshold.free}, then 
#' \code{f.criterion} is automatically set to \code{NULL}.
#' @param recall.levels a vector with the desired recall levels (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}) to compute the 
#' Precision at fixed Recall level (PXR). If \code{compute.performance=FALSE} the parameter \code{recall.levels} is automatically set to \code{NULL}.
#' @param compute.performance boolean value: should the flat and hierarchical performance (\code{AUPRC}, \code{AUROC}, \code{PXR}, 
#' \code{multilabel F-score}) be returned?    
#' \itemize{
#' \item \code{FALSE} (\code{def.}): performance are not computed and just the hierarchical scores matrix is returned;
#' \item \code{TRUE}: both performance and hierarchical scores matrix are returned;
#' }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension).
#' @param ann.file name of the file containing the label matrix of the examples (without rda extension).
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension).
#' @param flat.dir relative path where flat scores matrix is stored.
#' @param ann.dir relative path where annotation matrix is stored.
#' @param dag.dir relative path where graph is stored.
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored.
#' @param perf.dir relative path where the performance measures must be stored. If \code{compute.performance=FALSE} the functions 
#' automatically sets \code{perf.dir} to \code{NULL}.
#' @return Two \code{rda} files stored in the respective output directories:
#' \enumerate{
#'     \item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#'     for each example and for each considered class. It is stored in the \code{hierScore.dir} directory.
#'     \item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performance results:
#'     \enumerate{
#'         \item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'        \item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#'         \item PXR results computed though \code{precision.at.given.recall.levels.over.classes} (\code{\link{PXR}});
#'         \item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
#' }}
#' It is stored in the \code{perf.dir} directory.
#' @seealso \code{\link{GPAV}}
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' tmpdir <- paste0(tempdir(),"/");
#' save(g, file=paste0(tmpdir,"graph.rda"));
#' save(L, file=paste0(tmpdir,"labels.rda"));
#' save(S, file=paste0(tmpdir,"scores.rda"));
#' dag.dir <- flat.dir <- ann.dir <- tmpdir;
#' hierScore.dir <- perf.dir <- tmpdir;
#' recall.levels <- seq(from=0.25, to=1, by=0.25);
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' Do.GPAV(norm=FALSE, norm.type= "MaxNorm", W=NULL, parallel=FALSE, ncores=1, folds=NULL, 
#' seed=23, n.round=3, f.criterion ="F", recall.levels=recall.levels, compute.performance=TRUE, 
#' flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, 
#' dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.GPAV <- function(norm=TRUE, norm.type=NULL, W=NULL, parallel=FALSE, ncores=1, folds=5, seed=23, n.round=3, f.criterion ="F", 
    recall.levels=seq(from=0.1, to=1, by=0.1), compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
    dag.file=dag.file, flat.dir=flat.dir,ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir){
    
    ## Setting Check
    if(norm==FALSE && is.null(norm.type))
        stop("GPAV: If norm is set to FALSE, you need to specify a normalization method among those available", call.=FALSE);
    if(norm==TRUE && !is.null(norm.type))
        warning("GPAV: If norm is set to TRUE, the input flat matrix is already normalized.", paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
    if(parallel==TRUE && ncores<2)
        warning("GPAV: set ncores greater than 2 to exploit the GPAV parallel version", call.=FALSE);
    if(parallel==FALSE && ncores>=2)
        warning("GPAV: no GPAV parallel version is running, but ncores is higher or equal to 2.", " Set 'ncores' to 1 to run the sequential version or set 'parallel' to TRUE to run the parallel version", call.=FALSE);
    if(f.criterion!="F" && f.criterion!="avF" && compute.performance==TRUE)
        stop("GPAV: value of parameter 'f.criterion' misspelled", call.=FALSE);  
    if(compute.performance==FALSE && (!is.null(recall.levels) || !is.null(perf.dir) || !is.null(seed) || !is.null(n.round) || !is.null(f.criterion))){
        perf.dir <- NULL;
        recall.levels <- NULL;
        seed <- NULL;
        n.round <- NULL;
        f.criterion <- NULL;
    }
            
    ## loading dag
    dag.path <- paste0(dag.dir, dag.file,".rda");
    g <- get(load(dag.path));
    root <- root.node(g);

    ## loading annotation matrix
    ann.path <- paste0(ann.dir, ann.file,".rda");
    ann <- get(load(ann.path));
    ## removing root node from annotation table if it exists
    if(root %in% colnames(ann))
        ann <- ann[,-which(colnames(ann)==root)];

    ## loading flat matrix
    flat.path <- paste0(flat.dir, flat.file,".rda");
    if(norm){
        S <- get(load(flat.path));
        if(root %in% colnames(S)){
            root.scores <- S[,which(colnames(S)==root)];  ## needed to compute GPAV 
            S <- S[,-which(colnames(S)==root)];
        }    
    }else{
        S <- get(load(flat.path));
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "NORMALIZATION: DONE\n");
        if(root %in% colnames(S)){
            root.scores <- S[,which(colnames(S)==root)];  ## needed to compute GPAV 
            S <- S[,-which(colnames(S)==root)];
        }
    }

    ## check if |flat matrix classes| = |annotation matrix classes| 
    ## if not the classes of annotation matrix are shrunk to those of flat matrix
    class.check <- ncol(S)!=ncol(ann);
    if(class.check){
        ann <- ann[,colnames(S)];
        nd <- c(root,colnames(S));
        g <- do.subgraph(nd, g, edgemode="directed");
    }

    ## Compute FLAT PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        PRC.flat <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.flat <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.flat <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.flat <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("FLAT PERFORMANCE: DONE\n");
    }
    ## Hierarchical Correction 
    ## before running GPAV we need to re-add the root node scores
    if(!(root %in% colnames(S))){
        if(exists("root.scores")){
            S <- cbind(root.scores, S);
            colnames(S)[1] <- root;
        }else{
            max.score <- max(S);
            z <- rep(max.score, nrow(S));
            S <- cbind(z,S);
            colnames(S)[1] <- root;
        }
    }
    if(parallel){
        S <- GPAV.parallel(S, W=W, g, ncores=ncores);
    }else{
        S <- GPAV.over.examples(S, W=W, g);
    }
    S.hier <- S; # store hierarchical scores matrix with root node
    cat("HIERARCHICAL CORRECTION: DONE\n");

    ## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        ## remove root node before computing performance
        S <- S[,-which(colnames(S)==root)];
        PRC.hier <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.hier <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.hier <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.hier <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("HIERARCHICAL PERFORMANCE: DONE\n");
    }
    ## Storing Results
    rm(S);
    if(norm){
        save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        }
    }else{
        save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", norm.type, ".", flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        }
    }
}

#' @title GPAV holdout
#' @description High level function to correct the computed scores in a hierarchy according to the \code{GPAV} algorithm applying a
#' classical holdout procedure.
#' @details The function checks if the number of classes between the flat scores matrix and the annotations matrix mismatched.
#' If so, the number of terms of the annotations matrix is shrunk to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well. N.B.: it is supposed that all the nodes of the subgraph are accessible from the root.
#' @param norm boolean value:
#' \itemize{
#' \item \code{TRUE} (\code{def.}): the flat scores matrix has been already normalized in according to a normalization method;
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm} for which normalization can be applied;
#' }
#' @param norm.type can be one of the following three values:
#'  \enumerate{
#'  \item \code{NULL} (\code{def.}): set \code{norm.type} to \code{NULL} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used;
#'  }
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (\code{def. W=NULL}), \code{W} is a unitary
#' vector of the same length of the columns' number of the flat scores matrix (root node included).
#' @param parallel boolean value:
#' \itemize{
#'    \item \code{TRUE}: execute the parallel implementation of \code{GPAV} (\code{\link{GPAV.parallel}});
#'    \item \code{FALSE} (\code{def.}): execute the sequential implementation of \code{GPAV} (\code{\link{GPAV.over.examples}});
#' }
#' @param ncores number of cores to use for parallel execution (\code{def. 8}). Set the parameter \code{ncores} to \code{1} if the
#' parameter \code{parallel} is set to \code{FALSE}, otherwise set the desired number of cores.
#' @param folds number of folds of the cross validation on which computing the performance metrics averaged across folds (\code{def. 5}).
#' If \code{folds=NULL}, the performance metrics are computed one-shot, otherwise the performance metrics are averaged across folds.
#' If \code{compute.performance} is set to \code{FALSE}, \code{folds} is automatically set to \code{NULL}.
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. The parameter \code{seed} controls both the parameter \code{kk} and the parameter \code{folds}.
#' If \code{compute.performance} is set to \code{FALSE} and \code{bottomup} is set to \code{threshold.free}, then 
#' \code{seed} is automatically set to \code{NULL}.
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure.
#' If \code{compute.performance} is set to \code{FALSE} and \code{bottomup} is set to \code{threshold.free}, then 
#' \code{n.round} is automatically set to \code{NULL}.
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (def.): corresponds to the harmonic mean between the average precision and recall;
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples;
#' }
#' If \code{compute.performance} is set to \code{FALSE} and \code{bottomup} is set to \code{threshold.free}, then 
#' \code{f.criterion} is automatically set to \code{NULL}.
#' @param recall.levels a vector with the desired recall levels (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}) to compute the 
#' Precision at fixed Recall level (PXR). If \code{compute.performance=FALSE} the parameter \code{recall.levels} is automatically set to \code{NULL}.
#' @param compute.performance boolean value: should the flat and hierarchical performance (\code{AUPRC}, \code{AUROC}, \code{PXR}, 
#' \code{multilabel F-score}) be returned?    
#' \itemize{
#' \item \code{FALSE} (\code{def.}): performance are not computed and just the hierarchical scores matrix is returned;
#' \item \code{TRUE}: both performance and hierarchical scores matrix are returned;
#' }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension).
#' @param ann.file name of the file containing the label matrix of the examples (without rda extension).
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension).
#' @param ind.test.set name of the file containing a vector of integer numbers corresponding to the indices of the elements (rows) of scores
#' matrix to be used in the    test set.
#' @param ind.dir relative path to folder where \code{ind.test.set} is stored.
#' @param flat.dir relative path where flat scores matrix is stored.
#' @param ann.dir relative path where annotation matrix is stored.
#' @param dag.dir relative path where graph is stored.
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored.
#' @param perf.dir relative path where the performance measures must be stored. If \code{compute.performance=FALSE}, 
#' the parameter \code{perf.dir} is automatically set to \code{NULL}.
#' @return Two \code{rda} files stored in the respective output directories:
#' \enumerate{
#'     \item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#'     for each example and for each considered class. It is stored in the \code{hierScore.dir} directory;
#'     \item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performance results:
#'     \enumerate{
#'         \item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'         \item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#'         \item PXR results computed though \code{precision.at.given.recall.levels.over.classes} (\code{\link{PXR}});
#'         \item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
#' }}
#' It is stored in the \code{perf.dir} directory.
#' @seealso \code{\link{GPAV}}
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' data(test.index);
#' tmpdir <- paste0(tempdir(),"/");
#' save(g, file=paste0(tmpdir,"graph.rda"));
#' save(L, file=paste0(tmpdir,"labels.rda"));
#' save(S, file=paste0(tmpdir,"scores.rda"));
#' save(test.index, file=paste0(tmpdir,"test.index.rda"));
#' ind.dir <- dag.dir <- flat.dir <- ann.dir <- tmpdir;
#' hierScore.dir <- perf.dir <- tmpdir;
#' ind.test.set <- "test.index";
#' recall.levels <- seq(from=0.25, to=1, by=0.25);
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' Do.GPAV.holdout(norm=FALSE, norm.type="MaxNorm", W=NULL, parallel=FALSE, ncores=1, 
#' n.round=3, f.criterion="F", folds=NULL, seed=23, recall.levels=recall.levels, 
#' compute.performance=TRUE, flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
#' ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, 
#' dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.GPAV.holdout <- function(norm=TRUE, norm.type=NULL, W=NULL, parallel=FALSE, ncores=1, folds=5, seed=23, 
    n.round=3, f.criterion ="F", recall.levels=seq(from=0.1, to=1, by=0.1), compute.performance=FALSE, 
    flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, 
    flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir){
    
    ## Setting Check
    if(norm==FALSE && is.null(norm.type))
        stop("GPAV: If norm is set to FALSE, you need to specify a normalization method among those available", call.=FALSE);
    if(norm==TRUE && !is.null(norm.type))
        warning("GPAV: If norm is set to TRUE, the input flat matrix is already normalized.", paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
    if(parallel==TRUE && ncores<2)
        warning("GPAV: set ncores greater than 2 to exploit the GPAV parallel version", call.=FALSE);
    if(parallel==FALSE && ncores>=2)
        warning("GPAV: no GPAV parallel version is running, but ncores is higher or equal to 2.", " Set 'ncores' to 1 to run the sequential version or set 'parallel' to TRUE to run the parallel version", call.=FALSE);
    if(f.criterion!="F" && f.criterion!="avF" && compute.performance==TRUE)
        stop("GPAV: value of parameter 'f.criterion' misspelled", call.=FALSE);  
    if(compute.performance==FALSE && (!is.null(recall.levels) || !is.null(perf.dir) || !is.null(seed) || !is.null(n.round) || !is.null(f.criterion))){
        perf.dir <- NULL;
        recall.levels <- NULL;
        seed <- NULL;
        n.round <- NULL;
        f.criterion <- NULL;
    }

    ## Loading Data
    ## loading examples indices of the test set
    ind.set <- paste0(ind.dir, ind.test.set, ".rda");
    ind.test <- get(load(ind.set));

    ## loading dag
    dag.path <- paste0(dag.dir, dag.file,".rda");
    g <- get(load(dag.path));
    root <- root.node(g);

    ## loading annotation matrix
    ann.path <- paste0(ann.dir, ann.file,".rda");
    ann <- get(load(ann.path));
    ## removing root node from annotation table if it exists
    if(root %in% colnames(ann))
        ann <- ann[,-which(colnames(ann)==root)];

    ## loading flat matrix
    flat.path <- paste0(flat.dir, flat.file,".rda");
    if(norm){
        S <- get(load(flat.path));
        if(root %in% colnames(S)){
            root.scores <- S[,which(colnames(S)==root)];  ## needed to compute GPAV 
            S <- S[,-which(colnames(S)==root)];
        }
    }else{
        S <- get(load(flat.path));
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "NORMALIZATION: DONE", "\n");
        if(root %in% colnames(S)){
            root.scores <- S[,which(colnames(S)==root)];  ## needed to compute GPAV 
            S <- S[,-which(colnames(S)==root)];
        }
    }

    ## check if |flat matrix classes| = |annotation matrix classes| 
    ## if not the classes of annotation matrix are shrunk to those of flat matrix
    class.check <- ncol(S)!=ncol(ann);
    if(class.check){
        ann <- ann[,colnames(S)];
        nd <- c(root, colnames(S));
        g <- do.subgraph(nd, g, edgemode="directed");
    }

    ## shrinking scores flat matrix and annotation table to test test
    S <- S[ind.test,];    
    ann <- ann[ind.test,];

    ## Compute FLAT PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        PRC.flat <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.flat <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.flat <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.flat <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("FLAT PERFORMANCE: DONE", "\n");
    }
    ## Hierarchical Correction 
    ## before running GPAV we need to re-add the root node
    if(!(root %in% colnames(S))){
        if(exists("root.scores")){
            S <- cbind(root.scores[ind.test], S);
            colnames(S)[1] <- root;
        }else{
            max.score <- max(S);
            z <- rep(max.score, nrow(S));
            S <- cbind(z,S);
            colnames(S)[1] <- root;
        }
    }
    if(parallel){
        S <- GPAV.parallel(S, W=W, g, ncores=ncores);
    }else{
        S <- GPAV.over.examples(S, W=W, g);
    }
    ## before computing performance we need to remove the root node
    S <- S[,-which(colnames(S)==root)];
    cat("HIERARCHICAL CORRECTION: DONE", "\n");

    ## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        PRC.hier <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.hier <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.hier <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.hier <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("HIERARCHICAL PERFORMANCE: DONE", "\n");
    }
    ## Storing Results 
    S.hier <- S;
    rm(S);    
    if(norm){
        save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        }
    }else{
        save(S.hier, file=paste0(hierScore.dir, norm.type, ".", flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", norm.type, ".", flat.file, ".hierScores.GPAV.rda"), compress=TRUE);
        }
    }
}
