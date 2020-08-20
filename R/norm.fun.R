##########################
## Normalization Method ##
##########################

#' @title Max normalization
#' @description Function to normalize the scores of the flat scores matrix \code{S} per class.
#' @details the score of each class is normalized by dividing the score values for the maximum score of that class.
#' @param S matrix with the raw non normalized scores. Rows are examples and columns are classes.
#' @return A score matrix with the same dimensions of \code{S}, but with scores normalized.
#' @export
#' @examples
#' data(scores);
#' maxnorm <- normalize.max(S);
normalize.max <- function(S){
    classes <- colnames(S);
    maximum <- apply(S,2,max);
    for(class in classes){
        if(maximum[class] != 0){
            S[,class] <- S[,class]/maximum[class];
        }
    }
    return(S);
}

#' @title Scores Normalization Function
#' @description Functions to normalize a flat scores matrix w.r.t. max normalization (\code{maxnorm}) or quantile normalization (\code{qnorm}) 
#' @details To apply the quantile normalization the \pkg{preprocessCore} package must be properly installed.
#' @param norm.type can be one of the following two values:
#' \itemize{
#' \item maxnorm (\code{def.}): each score is divided w.r.t. the max of each class;
#' \item qnorm: a quantile normalization is applied. Package preprocessCore is used;
#' }
#' @param S a named flat scores matrix with examples on rows and classes on columns.
#' @return the matrix of the scores flat normalized w.r.t. \code{maxnorm} or \code{qnorm}.
#' @export
#' @examples
#' data(scores);
#' norm.types <- c("maxnorm","qnorm");
#' for(norm.type in norm.types){
#'     scores.normalization(norm.type=norm.type, S=S)
#' }
scores.normalization <- function(norm.type="maxnorm", S){
    if(norm.type=="maxnorm"){
        ## Max Normalization
        S <- normalize.max(S);        
    }else if(norm.type=="qnorm"){
        ## Quantile Normalization 
        ## NOTE: normalize.quantiles function returns a unnamed matrix. colnames are essential for hier.corr..
        S.norm <- normalize.quantiles(S);
        dimnames(S.norm) <- list(rownames(S),colnames(S));
        S <- S.norm;
        rm(S.norm);    
    }else{
        stop("scores.normalization: the chosen normalization method is not among those available or it has been misspelled");
    }
    return(S);
}
