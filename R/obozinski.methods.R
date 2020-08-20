###################################
##  Obozinski Heuristic Methods  ##
###################################

#' @name Heuristic-Methods
#' @aliases heuristic.max
#' @aliases heuristic.and
#' @aliases heuristic.or
#' @title Obozinski Heuristic Methods 
#' @description Implementation of the Heuristic Methods \code{MAX}, \code{AND}, \code{OR} (\cite{Obozinski et al., Genome Biology, 2008, 
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-s1-s6}{doi:10.1186/gb-2008-9-s1-s6}}).
#' @details Heuristic Methods:
#' \enumerate{
#'  \item \bold{MAX}: reports the largest logistic regression (LR) value of self and all descendants: \eqn{p_i = max_{j \in descendants(i)} \hat{p_j}};
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
        stop("heuristic.max: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);

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
        stop("heuristic.and: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);

    S.hier <- S;
    anc <- build.ancestors(g);
    for(i in 1:length(anc)){
      if(length(anc[[i]]) > 1){  # ancestors of i include also i
            curr.nd <- S[,names(anc[i])];
            forefathers <- as.matrix(S[,anc[[i]]]);
            curr.nd <- apply(forefathers, 1, prod);
            S.hier[,names(anc[i])] <- curr.nd;  
        }   
    }
    rm(S); gc();
    return(S.hier);
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
        stop("heuristic.or: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);

    S.hier <- S;
    desc <- build.descendants(g);
    for(i in 1:length(desc)){
      if(length(desc[[i]]) > 1){  # descendants of i include also i
            comp.progeny <- 1 - as.matrix(S[,desc[[i]]]);
            curr.nd <- apply(comp.progeny, 1, prod);        
            S.hier[,names(desc[i])] <- 1 - curr.nd;
        }
    }
    rm(S); gc();
    return(S.hier);
}

#' @title Call Heuristic Methods
#' @seealso \code{\link{Heuristic-Methods}}
#' @description Function to compute the hierarchical heuristic methods MAX, AND, OR (Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., Genome Biology, 2008}).
#' @param S a named flat scores matrix with examples on rows and classes on columns.
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param heuristic can be one of the following three values:
#' \enumerate{
#'  \item "MAX": run the method \code{heuristic.max};
#'  \item "AND": run the method \code{heuristic.and};
#'  \item "OR": run the method \code{heuristic.or};
#' }
#' @param norm boolean value. Should the flat score matrix be normalized? By default \code{norm=FALSE}. If \code{norm=TRUE} the matrix \code{S} is normalized according to \code{norm.type}.
#' @param norm.type can be one of the following values: 
#'  \enumerate{
#'  \item \code{NULL} (def.): none normalization is applied (\code{norm=FALSE})
#'  \item \code{maxnorm}: each score is divided for the maximum value of each class;
#'  \item \code{qnorm}: quantile normalization. \pkg{preprocessCore} package is used; 
#'  }
#' @return a matrix with the scores of the classes corrected according to the chosen heuristic algorithm.
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' S.and <- heuristic.methods(S, g, heuristic="and", norm=TRUE, norm.type="maxnorm");
heuristic.methods <- function(S, g, heuristic="and", norm=FALSE, norm.type=NULL){
    ## check
    if(heuristic!="max" && heuristic!="and" && heuristic!="or")
        stop("heuristic.methods: the chosen heuristic method is not among those available or it has been misspelled", call.=FALSE);
    if(norm==TRUE && is.null(norm.type))
        stop("heuristic.methods: choose a normalization methods among those available", call.=FALSE);
    if(norm==FALSE && !is.null(norm.type))
        warning("heuristic.methods: ", paste0("set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
   
    ## normalization
    if(norm){
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "normalization: done", "\n");
    }

    ## root node 
    root <- root.node(g);

    ## Obozinski's hierarchical heuristic methods 
    if(heuristic.fun=="and")
        S <- heuristic.and(S, g, root);
    if(heuristic.fun=="max")
        S <- heuristic.max(S, g, root);
    if(heuristic.fun=="or")
        S <- heuristic.or(S, g, root);  
    cat("Obozinski's heuristic", toupper(heuristic), "correction: done", "\n");
    return(S);
}

#' @title Holdout Heuristic Methods
#' @description Function to compute the hierarchical heuristic methods MAX, AND, OR (Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., Genome Biology, 2008}) applying a classical holdout procedure.
#' @param S a named flat scores matrix with examples on rows and classes on columns.
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @param testIndex a vector of integer numbers corresponding to the indexes of the elements (rows) of the scores matrix \code{S} to be used in the test set.
#' @param heuristic can be one of the following three values:
#' \enumerate{
#'  \item "max": run the method \code{heuristic.max};
#'  \item "and": run the method \code{heuristic.and};
#'  \item "or": run the method \code{heuristic.or};
#' }
#' @param norm boolean value. Should the flat score matrix be normalized? By default \code{norm=FALSE}. If \code{norm=TRUE} the matrix \code{S} is normalized according to \code{norm.type}.
#' @param norm.type can be one of the following values: 
#'  \enumerate{
#'  \item \code{NULL} (def.): none normalization is applied (\code{norm=FALSE})
#'  \item \code{maxnorm}: each score is divided for the maximum value of each class;
#'  \item \code{qnorm}: quantile normalization. \pkg{preprocessCore} package is used; 
#'  }
#' @return a matrix with the scores of the classes corrected according to the chosen heuristic algorithm. Rows of the matrix are shrunk to \code{testIndex}.
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(test.index);
#' S.and <- heuristic.holdout(S, g, testIndex=test.index, heuristic="and", norm=FALSE, norm.type=NULL);
heuristic.holdout <- function(S, g, testIndex, heuristic="and", norm=FALSE, norm.type=NULL){
    ## check
    if(heuristic!="max" && heuristic!="and" && heuristic!="or")
        stop("heuristic.methods: the chosen heuristic method is not among those available or it has been misspelled", call.=FALSE);
    if(norm==TRUE && is.null(norm.type))
        stop("heuristic.methods: choose a normalization methods among those available", call.=FALSE);
    if(norm==FALSE && !is.null(norm.type))
        warning("heuristic.methods: ", paste0("set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);

    ## normalization
    if(norm){
        S <- scores.normalization(norm.type=norm.type, S);
        cat(norm.type, "normalization: done", "\n");
    }

    ## root node 
    root <- root.node(g);

    ## shrinking scores matrix to test test
    S <- S[testIndex,];
   
    ## Obozinski's hierarchical heuristic methods 
    if(heuristic=="and")
        S <- heuristic.and(S, g, root);
    if(heuristic=="max")
        S <- heuristic.max(S, g, root);
    if(heuristic=="or")
        S <- heuristic.or(S, g, root);
    cat("Obozinski's heuristic", toupper(heuristic), "correction: done", "\n");
    return(S);
}


