###################################
##  Obozinski Heuristic Methods  ##
###################################

#' @name Heuristic-Methods
#' @aliases heuristic.max
#' @aliases heuristic.and
#' @aliases heuristic.or
#' @title Obozinski Heuristic Methods 
#' @description Implementation of the Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., Genome Biology, 2008, 
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-s1-s6}{doi:10.1186/gb-2008-9-s1-s6}}).
#' @details Heuristic Methods:
#' \enumerate{
#'  \item \bold{MAX}: reports the largest logist regression (LR) value of self and all descendants: \eqn{p_i = max_{j \in descendants(i)} \hat{p_j}};
#'  \item \bold{AND}: reports the product of LR values of all ancestors and self. This is equivalent to computing the probability that all 
#' ancestral terms are "on" assuming that, conditional on the data, all predictions are independent: \eqn{p_i = \prod_{j \in ancestors(i)} \hat{p_j}};
#'  \item \bold{OR}: computes the probability that at least one of the descendant terms is "on" assuming again that, conditional on the data, 
#' all predictions are independent: \eqn{1 - p_i = \prod_{j \in descendants(i)} (1 - \hat{p_j})};
#' }
#' @param S a named flat scores matrix with examples on rows and classes on columns.
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param root name of the class that it is the top-level (root) of the hierarchy (\code{def:00}).
#' @return a matrix with the scores of the classes corrected according to the chosen heuristic algorithm.
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' root  <- root.node(g);
#' S.max <- heuristic.max(S,g,root);
#' S.and <- heuristic.and(S,g,root);
#' S.or  <- heuristic.or(S,g,root);
heuristic.max <- function(S, g, root="00"){
    if(!(root %in% colnames(S))) {
        max.score <- max(S);
        z <- rep(max.score,nrow(S));
        S <- cbind(z,S);
        colnames(S)[1] <- root;
    }
    ## check consistency between nodes of g and classes of S
    class.check <- ncol(S)!=numNodes(g);
    if(class.check)
        stop("HEURISTIC: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);

    desc <- build.descendants(g);
    for(i in 1:length(desc)){
        curr.nd <- S[,names(desc[i])];
        progeny <- as.matrix(S[,desc[[i]]]);
        curr.nd <- apply(progeny, 1, max);
        S[,names(desc[i])] <- curr.nd;
    }
    return(S);
}

#' @rdname Heuristic-Methods
#' @export 
heuristic.and <- function(S, g, root="00"){
    if(!(root %in% colnames(S))) {
        max.score <- max(S);
        z <- rep(max.score,nrow(S));
        S <- cbind(z,S);
        colnames(S)[1] <- root;
    }
    ## check consistency between nodes of g and classes of S
    class.check <- ncol(S)!=numNodes(g);
    if(class.check)
        stop("HEURISTIC: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);

    S.corr <- S;
    anc <- build.ancestors(g);
    for(i in 1:length(anc)){
      if(length(anc[[i]]) > 1){  # ancestors of i include also i
            curr.nd <- S[,names(anc[i])];
            forefathers <- as.matrix(S[,anc[[i]]]);
            curr.nd <- apply(forefathers, 1, prod);
            S.corr[,names(anc[i])] <- curr.nd;  
        }   
    }
    return(S.corr);
}

#' @rdname Heuristic-Methods
#' @export 
heuristic.or <- function(S, g, root="00"){
    if(!(root %in% colnames(S))) {
        max.score <- max(S);
        z <- rep(max.score,nrow(S));
        S <- cbind(z,S);
        colnames(S)[1] <- root;
    }
    ## check consistency between nodes of g and classes of S
    class.check <- ncol(S)!=numNodes(g);
    if(class.check)
        stop("HEURISTIC: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);

    S.corr <- S;
    desc <- build.descendants(g);
    for(i in 1:length(desc)){
      if(length(desc[[i]]) > 1){  # descendants of i include also i
            comp.progeny <- 1 - as.matrix(S[,desc[[i]]]);
            curr.nd <- apply(comp.progeny, 1, prod);        
            S.corr[,names(desc[i])] <- 1 - curr.nd;
        }
    }
    return(S.corr);
}

#' @title Do Heuristic Methods
#' @seealso \code{\link{Heuristic-Methods}}
#' @description High level function to compute the hierarchical heuristic methods MAX, AND, OR (Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., 
#' Genome Biology, 2008}).
#' @details The function checks if the number of classes between the flat scores matrix and the annotations matrix mismatched.
#' If so, the number of terms of the annotations matrix is shrunk to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well. N.B.: it is supposed that all the nodes of the subgraph are accessible from the root.
#' @details We excluded the predictions of the root node in computing all the performances, since it is a \emph{dummy} node added 
#' to the ontology for practical reasons (e.g. some graph-based software may require a single root node to work). However, the root node scores 
#' are stored in the hierarchical scores matrix.
#' @param heuristic.fun can be one of the following three values:
#' \enumerate{
#'  \item "MAX": run the heuristic method MAX;
#'  \item "AND": run the heuristic method AND;
#'  \item "OR": run the heuristic method OR;
#' }
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method; 
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied;
#' }
#' @param norm.type can be one of the following three values: 
#'  \enumerate{
#'  \item \code{NULL} (def.): set \code{norm.type} to \code{NULL} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used; 
#'  }
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
#'  \item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#'  for each example and for each considered class. It is stored in the \code{hierScore.dir} directory;
#'  \item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performace results:
#'  \enumerate{
#'      \item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'      \item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#'      \item PXR results computed though \code{precision.at.given.recall.levels.over.classes} (\code{\link{PXR}});
#'      \item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
#' }}
#' It is stored in the \code{perf.dir} directory.
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
#' Do.heuristic.methods(heuristic.fun="AND", norm=FALSE, norm.type="MaxNorm",
#' folds=NULL, seed=23, n.round=3, f.criterion ="F", recall.levels=recall.levels, 
#' compute.performance=TRUE, flat.file=flat.file, ann.file=ann.file, 
#' dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.heuristic.methods <- function(heuristic.fun="AND", norm=TRUE, norm.type=NULL, folds=5, seed=23, 
    n.round=3, f.criterion ="F", recall.levels=seq(from=0.1, to=1, by=0.1), compute.performance=FALSE, 
    flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, 
    dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir){

    ## Setting Check
    if(heuristic.fun!="MAX" && heuristic.fun!="AND" && heuristic.fun!="OR")
        stop("HEURISTIC: the chosen heuristic method is not among those available or it has been misspelled", call.=FALSE);
    if(norm==FALSE && is.null(norm.type))
        stop("HEURISTIC: If norm is set to FALSE, you need also to specify a normalization method among those available", call.=FALSE);
    if(norm==TRUE && !is.null(norm.type))
        warning("HEURISTIC: If norm is set to TRUE, the input flat matrix is already normalized.", paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
    if(f.criterion!="F" && f.criterion!="avF" && compute.performance==TRUE)
        stop("HEURISTIC: value of parameter 'f.criterion' misspelled", call.=FALSE);  
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
    if(root %in% colnames(ann))
        ann <- ann[,-which(colnames(ann)==root)];

    ## loading flat matrix
    flat.path <- paste0(flat.dir, flat.file,".rda");
    if(norm){
        S <- get(load(flat.path));
        if(root %in% colnames(S)){
            classes.names <- colnames(S)
            root.scores <- S[,which(colnames(S)==root)];
            S <- S[,-which(colnames(S)==root)];
        }
    }else{
        S <- get(load(flat.path));
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "NORMALIZATION: DONE", "\n");
        if(root %in% colnames(S)){
            classes.names <- colnames(S)
            root.scores <- S[,which(colnames(S)==root)];
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

    ## Compute FLAT PRC, AUC, PXR (average and per class) FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        PRC.flat <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.flat <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.flat <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.flat <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("FLAT PERFORMANCE: DONE", "\n");
    }

    ## Obozinski's Hierarchical Heuristic Methods 
    if(exists("root.scores")){
        S <- cbind(root.scores, S);
        colnames(S)[1] <- root;
        S <- S[,classes.names];
    }
    if(heuristic.fun=="AND")
        S <- heuristic.and(S, g, root);
    if(heuristic.fun=="MAX")
        S <- heuristic.max(S, g, root);
    if(heuristic.fun=="OR")
        S <- heuristic.or(S, g, root);  
    S.hier <- S;  
    cat("HIERARCHICAL CORRECTION: DONE", "\n");

    ## Compute HIER PRC, AUC, PXR (average and per class) FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        if(root %in% colnames(S))
            S <- S[,-which(colnames(S)==root)];
        PRC.hier <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.hier <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.hier <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.hier <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("HIERARCHICAL PERFORMANCE: DONE", "\n");
    }   

    ## Storing Results #########
    rm(S); 
    if(heuristic.fun=="AND")
        heuristic.name <- "AND";
    if(heuristic.fun=="MAX")
        heuristic.name <- "MAX";
    if(heuristic.fun=="OR")
        heuristic.name <- "OR";
    if(norm){
        save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores",heuristic.name,".rda"), compress=TRUE);
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
        }
    }else{
        save(S.hier, file=paste0(hierScore.dir, norm.type, ".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);   
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", norm.type, ".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);   
        }
    }
}

#' @title Do Heuristic Methods holdout
#' @description High level function to compute the hierarchical heuristic methods MAX, AND, OR (Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., 
#' Genome Biology, 2008}) applying a classical holdout procedure.
#' @details The function checks if the number of classes between the flat scores matrix and the annotations matrix mismatched.
#' If so, the number of terms of the annotations matrix is shrunk to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well. N.B.: it is supposed that all the nodes of the subgraph are accessible from the root.
#' @details We excluded the predictions of the root node in computing all the performances, since it is a \emph{dummy} node added 
#' to the ontology for practical reasons (e.g. some graph-based software may require a single root node to work). However, the root node scores 
#' are stored in the hierarchical scores matrix.
#' @param heuristic.fun can be one of the following three values:
#' \enumerate{
#'  \item "MAX": run the heuristic method MAX;
#'  \item "AND": run the heuristic method AND;
#'  \item "OR": run the heuristic method OR;
#' }
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method; 
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm} for which normalization can be applied;
#' }
#' @param norm.type can be one of the following three values:
#'  \enumerate{
#'  \item \code{NONE} (def.): set \code{norm.type} to \code{NONE} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used;
#'  }
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
#' Precision at fixed Recall level (PXR). If \code{compute.performance=FALSE} the parameter \code{recall.levels} is automatically set to \code{NULL}
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
#' matrix to be used in the test set.
#' @param ind.dir relative path to folder where \code{ind.test.set} is stored.
#' @param flat.dir relative path where flat scores matrix is stored.
#' @param ann.dir relative path where annotation matrix is stored.
#' @param dag.dir relative path where graph is stored.
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored.
#' @param perf.dir relative path where the performance measures must be stored. If \code{compute.performance=FALSE}, 
#' the parameter \code{perf.dir} is automatically set to \code{NULL}.
#' @return Two \code{rda} files stored in the respective output directories:
#' \enumerate{
#'  \item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#'  for each example and for each considered class. It is stored in the \code{hierScore.dir} directory;
#'  \item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performace results:
#'  \enumerate{
#'      \item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'      \item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#'      \item PXR results computed though \code{precision.at.given.recall.levels.over.classes} (\code{\link{PXR}});
#'      \item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
#' }}
#' It is stored in the \code{perf.dir} directory.
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
#' recall.levels <- seq(from=0.25, to=1, by=0.25);
#' ind.test.set <- "test.index";
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' Do.heuristic.methods.holdout(heuristic.fun="MAX", norm=FALSE, norm.type="MaxNorm", 
#' folds=NULL, seed=23, n.round=3, f.criterion ="F", recall.levels=recall.levels,
#' compute.performance=TRUE, flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
#' ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, 
#' dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.heuristic.methods.holdout <- function(heuristic.fun="AND", norm=TRUE, norm.type=NULL, folds=5, seed=23, 
    n.round=3, f.criterion ="F", recall.levels=seq(from=0.1, to=1, by=0.1), compute.performance=FALSE, 
    flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, 
    flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir){

    ## Setting Check
    if(heuristic.fun!="MAX" && heuristic.fun!="AND" && heuristic.fun!="OR")
        stop("HEURISTIC: the chosen heuristic method is not among those available or it has been misspelled", call.=FALSE);
    if(norm==FALSE && is.null(norm.type))
        stop("HEURISTIC: If norm is set to FALSE, you need also to specify a normalization method among those available", call.=FALSE);
    if(norm==TRUE && !is.null(norm.type))
        warning("HEURISTIC: If norm is set to TRUE, the input flat matrix is already normalized.", paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
    if(f.criterion!="F" && f.criterion!="avF" && compute.performance==TRUE)
        stop("HEURISTIC: value of parameter 'f.criterion' misspelled", call.=FALSE);  
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
    if(root %in% colnames(ann))
        ann <- ann[,-which(colnames(ann)==root)];

     ## loading flat matrix
    flat.path <- paste0(flat.dir, flat.file,".rda");
    if(norm){
        S <- get(load(flat.path));
        if(root %in% colnames(S)){
            classes.names <- colnames(S)
            root.scores <- S[,which(colnames(S)==root)];
            S <- S[,-which(colnames(S)==root)];
        }
    }else{
        S <- get(load(flat.path));
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "NORMALIZATION: DONE", "\n");
        if(root %in% colnames(S)){
            classes.names <- colnames(S)
            root.scores <- S[,which(colnames(S)==root)];
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

    ## Compute FLAT PRC, AUC, PXR (average and per class) FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        PRC.flat <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.flat <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.flat <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.flat <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("FLAT PERFORMANCE: DONE", "\n");
    }
    ## Obozinski's Hierarchical Heuristic Methods 
    if(exists("root.scores")){
        S <- cbind(root.scores[ind.test], S);
        colnames(S)[1] <- root;
        S <- S[,classes.names];
    }
    if(heuristic.fun=="AND")
        S <- heuristic.and(S, g, root);
    if(heuristic.fun=="MAX")
        S <- heuristic.max(S, g, root);
    if(heuristic.fun=="OR")
        S <- heuristic.or(S, g, root);
    S.hier <- S; 
    cat("HIERARCHICAL CORRECTION: DONE", "\n");

    ## Compute HIER PRC, AUC, PXR (average and per class) FMM (average and per-example) one-shoot or cross-validated 
    if(compute.performance){
        if(root %in% colnames(S))
            S <- S[,-which(colnames(S)==root)];
        PRC.hier <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
        AUC.hier <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
        PXR.hier <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
        FMM.hier <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE, folds=folds, seed=seed);
        cat("HIERARCHICAL PERFORMANCE: DONE", "\n");
    }
    
    ## Storing Results 
    rm(S);
    if(heuristic.fun=="AND")
        heuristic.name <- "AND";
    if(heuristic.fun=="MAX")
        heuristic.name <- "MAX";
    if(heuristic.fun=="OR")
        heuristic.name <- "OR";

    if(norm){
        save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores",heuristic.name,".rda"), compress=TRUE);
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
        }
    }else{
        save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);    
        if(compute.performance){
            save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, file=paste0(perf.dir, "PerfMeas.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);    
        }
    }
}
