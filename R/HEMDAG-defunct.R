#' @name HEMDAG-defunct
#' @title Defunct functions in package \pkg{HEMDAG}.
#' @description The functions listed below are defunct. The alternative function is supplied.
#' Help page for defunct functions is available by typing \code{help("HEMDAG-defunct")}

##*********##
## DESCENS ##
##*********##
#' @param S a named flat scores matrix with examples on rows and classes on columns
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is on the top-level of the hierarchy (\code{def. root="00"})
#' @param t threshold for the choice of the positive descendants (\code{def. t=0.5}); whereas in the \code{descens.tau} variant 
#' the parameter \code{t} balances the contribution between the positives children of a node \eqn{i} and that of its
#' positives descendants excluding its positives children
#' @param w weight to balance between the contribution of the node \eqn{i} and that of its positive descendants
#' @return a named matrix with the scores of the classes corrected according to the chosen hierarchical variant.

#' @templateVar old descens.threshold
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export 
descens.threshold <- function(S, g, root="00", t=0.5){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > t;    # positive descendants selection
					desc.pos <- desc[j,][desc.set];
					parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
				}
				S[,names(desc.bup[i])] <- parent;
			}
		}
		# top-down visit
		par.tod <- get.parents.top.down(g,levels,root)
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @templateVar old descens.threshold.free
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export 
descens.threshold.free <- function(S, g, root="00"){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > parent[j];	# positive descendants selection
					desc.pos <- desc[j,][desc.set];
					parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
				}
				S[,names(desc.bup[i])] <- parent;
			}
		}
		# top-down visit
		par.tod <- get.parents.top.down(g,levels,root)
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @templateVar old descens.weighted.threshold.free
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export 
descens.weighted.threshold.free <- function(S, g, root="00", w=0.5){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > parent[j];
					desc.pos <- desc[j,][desc.set];
					if(length(desc.pos)!=0){
						parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
					}
				}
				S[,names(desc.bup[i])] <- parent;
			}
		}
		# top-down visit
		par.tod <- get.parents.top.down(g,levels,root)
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @templateVar old descens.weighted.threshold
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export 
descens.weighted.threshold <- function(S, g, root="00", t=0.5, w=0.5){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > t;
					desc.pos <- desc[j,][desc.set];
					if(length(desc.pos)!=0){
						parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
					}
				}
				S[,names(desc.bup[i])] <- parent;
			}
		}
		# top-down visit
		par.tod <- get.parents.top.down(g,levels,root)
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @templateVar old descens.tau
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export 
descens.tau <- function(S, g, root="00", t=0.5){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
	for(i in 1:length(desc.bup)){
		if(length(desc.bup[[i]])!=1){ 		# in the desc list is included also the node itself
			node.curr <- nodes[i];
			parent <- S[,names(desc.bup[i])];
			tmp <- setdiff(desc.bup[[i]], node.curr);		# descendants
			delta <- setdiff(tmp, chd.bup[[i]]);  			# descendants without children 
			children <- as.matrix(S[,chd.bup[[i]]]);		# genes considering children node 
			desc <-  as.matrix(S[,delta]);					# genes considering descendants nodes without children
			for(j in 1:length(parent)){
				desc.set <- desc[j,] > parent[j];			# positive descendants (without children) selection
				desc.pos <- desc[j,][desc.set];
				child.set <- children[j,] > parent[j];  	# positive children selection
				child.pos <- children[j,][child.set];
				parent[j] <- t * ((parent[j] + sum(child.pos))/(1+length(child.pos))) + (1- t) * ((parent[j] + sum(desc.pos))/(1+length(desc.pos)));
			}
			S[,names(desc.bup[i])] <- parent;
		}
	}
	# top-down visit
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

##*********##
## TPR-DAG ##
##*********##
#' @templateVar old tpr.threshold
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export
tpr.threshold <- function(S, g, root="00", t=0.5){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("TPR-DAG: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > t;    # positive children selection
				child.pos <- children[j,][child.set];
				parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat scores correction
			}
			S[,names(chd.bup[i])] <- parent;
		}
	}
	# top-down visit
	par.tod <- get.parents.top.down(g,levels,root)
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @templateVar old tpr.threshold.free
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export
tpr.threshold.free <- function(S, g, root="00"){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("TPR-DAG: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels	
	levels <- graph.levels(g,root);
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > parent[j]; # positive children selection
				child.pos <- children[j,][child.set];
				parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat score correction
			}
			S[,names(chd.bup[i])] <- parent;
		} 
	}
	# top-down visit
	par.tod <- get.parents.top.down(g,levels,root);
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
		parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
		# Note: the version with an apply and an ifelse statement is slower ...
		for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
				child[j] <- x;   # hierarchical correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @templateVar old tpr.weighted.threshold.free
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export
tpr.weighted.threshold.free <- function(S, g, root="00", w=0.5){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("TPR-DAG: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > parent[j];    # positive children selection
				child.pos <- children[j,][child.set];
				if(length(child.pos)!=0){
					parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score correction
				}
			}
			S[,names(chd.bup[i])] <- parent;
		}
	}
	# top-down visit
	par.tod <- get.parents.top.down(g,levels,root)
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @templateVar old tpr.weighted.threshold
#' @templateVar new TPR.DAG
#' @template template-defun_pkg
#' @export
tpr.weighted.threshold <- function(S, g, root="00", t=0.5, w=0.5){
	.Defunct("TPR.DAG");
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("TPR-DAG: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > t;    # positive children selection
				child.pos <- children[j,][child.set];
				if(length(child.pos)!=0){
					parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score prediction
				}
			}
			S[,names(chd.bup[i])] <- parent;
		}
	}
	# top-down visit grafo
	par.tod <- get.parents.top.down(g,levels,root)
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
	S <- S[,-which(colnames(S)==root)];
	return(S);
}
