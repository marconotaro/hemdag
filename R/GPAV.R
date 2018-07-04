##******##
## GPAV ##
##******##

#' @title Binary Square Adjacency Matrix
#' @description This function returns a binary square upper triangular matrix where rows and colums correspond to the nodes's name of the graph \code{g}. 
#' @details The nodes of the matrix are topologically sorted. Let's denote with \code{adj} our adjacency matrix. Then \code{adj} represents a partial 
#' order data set in which the class \code{j} dominates the class \code{i}. In other words, \code{adj[i,j]=1} means that \code{j} dominates \code{i}; 
#' \code{adj[i,j]=0} means that there is no edge between the class \code{i} and the class \code{j}. Moreover the nodes of \code{adj} are enumerated 
#' so that \code{adj[i,j]=1} implies \eqn{i < j}, i.e. \code{adj} is upper triangular.
#' @param g a graph of class \code{graphNELL} representing the hierarchy of the class.
#' @return an adjacency matrix which is square, logical and upper triangular 
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
#' @description Implementetion of GPAV (Generalized Pool-Adjacent Violators) algorithm 
#' (\cite{Burdavok et al., Journal of Computational Mathematics,  2006 -- \href{http://www.jstor.org/stable/43693336}{link}}) 
#' to correct the scores of the hierarchy according to the constraints that the score of a node cannot be greater than a score of its parents.
#' GPAV algorithm treats nodes sequentially (in accordance with the topological order of \code{adj}) and merges some adjacent nodes into sets
#' which are called blocks. The fitted values for the nodes of the same block are equal to the weighted average value of their response values.
#' @details Given the constraints adjacency matrix of the graph, a vector of scores \eqn{\hat{y} \in R^n} and a vector of strictly positive 
#' weights \eqn{w \in R^n}, the GPAV algorithm returns a vector \eqn{\bar{y}} which is as close as possible, in the least-squares sense, 
#' to the response vector \eqn{\hat{y}} and whose components are partially ordered in accordance with the constraints matrix \code{adj}. 
#' In other words, GPAV solves the following problem:
#' \deqn{
#'	\bar{y} = \left\{
#'	\begin{array}{l}
#'		\min \sum_{i \in N} (\hat{y}_i - \bar{y}_i )^2\\\\
#'		\forall i, \quad  j \in par(i) \Rightarrow  \bar{y}_j  \geq \bar{y}_i
#'	\end{array}
#' \right.
#'}
#' @param Y vector of scores relative to a single example. \code{Y} must be a numeric named vector, where names
#' correspond to classes' names, i.e. nodes of the graph \code{g} (root node included)
#' @param W vector of weight relative to a single example. Vector \code{Y} and \code{W} must have the same length
#' @param adj adjacency matrix of the graph which is sparse, logical and upper triangular. Number of colums of \code{adj} must be 
#' equal to the length of \code{Y} and \code{W}
#' @seealso \code{\link{adj.upper.tri}}
#' @return a list of 3 elements:
#' \itemize{
#'	\item \code{YFit}: a named vector with the scores of the classes corrected according to the GPAV algorithm. 
#'	\code{NOTE}: the classes of \code{YFit} are topologically sorted, that is are in the same order of those of \code{adj}.
#'	\item \code{blocks}: list of vectors, containing the partitioning of nodes (represented with an integer number) into blocks;
#'	\item \code{W}: vector of weights.
#' }
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' Y <- S[3,];
#' W <- rep(1,ncol(S));
#' adj <- adj.upper.tri(g);
#' S.GPAV <- GPAV(Y,W,adj);
GPAV <- function(Y, W, adj){
	nW <- length(W);
	nY <- length(Y);
	nadj <- ncol(adj);
	if(nY!=nadj){
		stop("GPAV: number of classes between 'Y' and 'adj' does not match.");
	}
	if(nW!=nadj){
		stop("GPAV: number of classes between 'W' and 'adj' does not match.");
	}

	## sort nodes in the same order of adj matrix (i.e. in a topologically ordered)
	Y <- Y[colnames(adj)];
	
	## just for cleannes: we play with index and not with name...
	Y <- unname(Y);
	W <- unname(W);

	## assign to each term an integer number, i.e. index's term
	N <- ncol(adj);
	corr <- 1:N;

	## initialize absorbing blocks;
	blocks <- vector(mode="list", length=N);
	for(j in 1:N){
		blocks[[j]] <- j;
	}

	## GPAV algorithm	
	for(j in 1:N){
		Items <- 1;
		while(length(Items)!=0){
			## index of j's children/child
			if(is.matrix(adj[,blocks[[j]]])){
				jMinus <- unname(which(adj[,blocks[[j]]]==1, arr.ind=TRUE)[,1]);
			}else{
				jMinus <- unname(which(adj[,blocks[[j]]]==1));
			}
			jMinus <- unique(corr[jMinus]); ## 'unique' reduces redundance due to block absorbing..
			Items <- which((Y[jMinus] > Y[j])==TRUE); ## constraints' violation: children > parent
			if(length(Items)!=0){
				LookingIn <- jMinus[Items];  
				Xmax <- max(Y[LookingIn]);
				i.xmax <- which(Xmax==Y[LookingIn]);
				i <- LookingIn[i.xmax][1];
				## NOTE: if we deal with diff equal Xmax, we take the 1st one; 
				## the remaining will be considered in the 2nd while loop..
				
				## Merging
				Y[j] <- (W[i]*Y[i] + W[j]*Y[j])/(W[j] + W[i]);
				W[j] <- W[j] + W[i];
				
				## Add block i into j
				blocks[[j]] <- append(unlist(blocks[i]), blocks[[j]]);
				blocks[[i]] <- 0;
				corr[blocks[[j]]] <- j;
				Y[i] <- 0;
			}
		}
	}
	YFit <- Y[corr];
	names(YFit) <- colnames(adj);
	names(W) <- colnames(adj);
	res <- list(YFit=YFit, blocks=blocks, W=W);
	return(res);
}

#' @title GPAV Over Examples
#' @description Function to compute GPAV across all the examples
#' @param S a named flat scores matrix with examples on rows and classes on columns (root node included)
#' @param W vector of weight relative to a single example. The length of \code{W} must be equal to the number of nodes of \code{g} or to the 
#' colums' number of \code{S} (including the root node)
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @return a named matrix with the scores of the classes corrected according to the GPAV algorithm
#' @seealso \code{\link{GPAV.parallel}}
#' @export 
#' @examples
#' data(graph);
#' data(scores);
#' W <- rep(1,ncol(S));
#' S.GPAV <- GPAV.over.examples(S,W,g);
GPAV.over.examples <- function(S, W, g){
	adj <- adj.upper.tri(g);
	M <- c();
	for(i in 1:nrow(S)){
		M <- rbind(M, GPAV(S[i,],W,adj)$YFit);
	}
	rownames(M) <- rownames(S);
	M <- M[,colnames(S)];
	S <- M;	rm(M);
	return(S);
}

#' @title GPAV Over Examples -- Parallel Implementation
#' @description Function to compute GPAV across all the examples (parallel implementation)
#' @param S a named flat scores matrix with examples on rows and classes on columns (root node included)
#' @param W vector of weight relative to a single example. The length of \code{W} must be equal to the number of nodes of \code{g} or to the 
#' colums' number of \code{S} (including the root node)
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param ncores number of cores to use for parallel execution (\code{def. 8})
#' @return a named matrix with the scores of the classes corrected according to the GPAV algorithm
#' @export 
#' @examples
#' data(graph);
#' data(scores);
#' W <- rep(1,ncol(S));
#' if (Sys.info()['sysname']!="Windows"){
#'	S.GPAV <- GPAV.parallel(S,W,g,ncores=2);
#' }
GPAV.parallel <- function(S, W, g, ncores=8){
	adj <- adj.upper.tri(g);
	if(ncores == 0){
		n.cores <- detectCores(); 
		if (n.cores > 3)
		ncores <- n.cores - 1;
	}
	registerDoParallel(cores=ncores);
	res.list <- foreach(i=1:nrow(S), .inorder=FALSE) %dopar% {
		res <- GPAV(S[i,],W,adj)$YFit;		
	}
	M <- do.call(rbind,res.list);
	rownames(M) <- rownames(S);
	M <- M[,colnames(S)];
	S <- M;	rm(M);
	return(S);
}
