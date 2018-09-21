##******************##
## TPR-DAG VARIANTS ##
##******************##
#' @name TPR-DAG-variants
#' @seealso \code{\link{GPAV}}, \code{\link{HTD-DAG}}
#' @title TPR-DAG Ensemble Variants
#' @description Function gathering the true-path-rule-based hierarchical learning ensemble algorithms and its variants. 
#' In their more general form the \code{TPR-DAG} algorithms adopt a two step learning strategy:
#' \enumerate{
#'	\item in the first step they compute a \emph{per-level bottom-up} visit from the leaves to the root to propagate positive predictions 
#'  across the hierarchy;
#'	\item in the second step they compute a \emph{per-level top-down} visit from the root to the leaves in order to assure the hierarchical 
#' 	consistency of the predictions
#' }
#' @details The \emph{vanilla} \code{TPR-DAG} adopts a per-level bottom-up traversal of the DAG to correct the flat predictions \eqn{\hat{y}_i}:
#' \deqn{
#' 	\bar{y}_i := \frac{1}{1 + |\phi_i|} (\hat{y}_i + \sum_{j \in \phi_i} \bar{y}_j)
#' }
#' where \eqn{\phi_i} are the positive children of \eqn{i}.
#' Different strategies to select the positive children \eqn{\phi_i} can be applied:
#' \enumerate{
#' 	\item \strong{Threshold-Free} strategy: the positive nodes are those children that can increment the score of the node \eqn{i}, that is those nodes 
#' 	that achieve a score higher than that of their parents:
#' 	\deqn{
#' 		\phi_i := \{ j \in child(i) | \bar{y}_j > \hat{y}_i \}
#' 	}
#' 	\item \strong{Threshold} strategy: the positive children are selected on the basis of a threshold that can ben selected in two different ways:
#' 	\enumerate{
#' 		\item for each node a constant threshold \eqn{\bar{t}} is a priori selected:
#'		\deqn{
#'			\phi_i := \{ j \in child(i) | \bar{y}_j > \bar{t} \}
#'		}
#' 		For instance if the predictions represent probabilities it could be meaningful to a priori select \eqn{\bar{t}=0.5}.
#' 		\item the threshold is selected to maximize some performance metric \eqn{\mathcal{M}} estimated on the training data, as for instance
#' 		the F-score or the AUPRC. In other words the threshold is selected to maximize some measure of accuracy of the predictions 
#' 		\eqn{\mathcal{M}(j,t)} on the training data for the class \eqn{j} with respect to the threshold \eqn{t}. 
#' 		The corresponding set of positives \eqn{\forall i \in V} is:
#' 		\deqn{
#' 			\phi_i := \{ j \in child(i) | \bar{y}_j > t_j^*,  t_j^* = \arg \max_{t} \mathcal{M}(j,t) \}
#' 		}
#' 		For instance \eqn{t_j^*} can be selected from a set of \eqn{t \in (0,1)} through internal cross-validation techniques.
#'	}
#' }
#' The weighted \code{TPR-DAG} version can be designed by adding a weight \eqn{w \in [0,1]} to balance between the 
#' contribution of the node \eqn{i} and that of its positive children \eqn{\phi}, through their convex combination:
#' \deqn{
#' 	\bar{y}_i := w \hat{y}_i + \frac{(1 - w)}{|\phi_i|} \sum_{j \in \phi_i} \bar{y}_j
#' }
#' If \eqn{w=1} no weight is attributed to the children and the \code{TPR-DAG} reduces to the \code{HTD-DAG} algorithm, since in this
#' way only the prediction for node \eqn{i} is used in the bottom-up step of the algorithm. If \eqn{w=0} only the predictors 
#' associated to the children nodes vote to predict node \eqn{i}. In the intermediate cases we attribute more importance to the predictor for the
#' node \eqn{i} or to its children depending on the values of \eqn{w}. \cr
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
#' \eqn{\delta_i} positive nodes. \cr
#'
#' Simply by replacing the \code{HTD} (\code{\link{HTD-DAG}}) top-down step with the \code{GPAV} approach (\code{\link{GPAV}}) we can desing the
#' \code{TPR-DAG} variant \code{ISO-TPR}. The most important feature of \code{ISO-TPR} is that it maintains the hierarchical constraints by
#' construction and selects the closest solution (in the sense of the least squared error) to the flat predictions that obeys to the true path rule.
#' Obviously, any aforementioned strategy for the selection of ``positive'' children or descendants can be applied before executing the \code{GPAV} correction.
#' @param S a named flat scores matrix with examples on rows and classes on columns
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is on the top-level of the hierarchy (\code{def. root="00"})
#' @param positive choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#' 	\item \code{children} (\code{def.}): for each node are considered its positive children;
#' 	\item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#' 	\item \code{threshold.free} (\code{def.}): positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#' 	\item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#' 	\item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#' 	\item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#' 	\item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy. 
#'	NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positive=descendants};
#' }
#' @param topdown strategy to make the scores hierarchy-consistent. It can be one of the following values:
#' \itemize{
#' 	\item \code{HTD} (\code{def.}): \code{HTD-DAG} strategy is applied (\code{\link{HTD-DAG}});
#' 	\item \code{GPAV}: \code{GPAV} strategy is applied (\code{\link{GPAV}}).
#' }
#' @param t threshold for the choice of positive nodes (\code{def. t=0}). Set \code{t} only for the variants that requiring 
#' a threshold for the selection of the positive nodes, otherwise set \code{t} to zero
#' @param w weight to balance between the contribution of the node \eqn{i} and that of its positive nodes. Set \code{w} only for the
#' \emph{weighted} variants, otherwise set \code{w} to zero
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (by \code{def. W=NULL}), \code{W} is a unitary 
#' vector of the same length of the columns' number of the flat scores matrix (root node included). Set \code{W} only if \code{topdown=GPAV}.
#' @param parallel boolean value:
#' \itemize{
#'	\item \code{TRUE}: execute the parallel implementation of GPAV (\code{\link{GPAV.parallel}});
#'	\item \code{FALSE} (def.): execute the sequential implementation of GPAV (\code{\link{GPAV.over.examples}}).
#' }
#' Use \code{parallel} if and only if \code{topdown=GPAV}; otherwise set \code{parallel=FALSE}.
#' @param ncores number of cores to use for parallel execution (\code{def. 8}). Set \code{ncores=1} if \code{parallel=FALSE}, 
#' otherwise set \code{ncores} to the desired number of cores.
#' Use \code{ncores} if and only if \code{topdown=GPAV}; otherwise set \code{parallel=1}.
#' @return a named matrix with the scores of the classes corrected according to the chosen algorithm
#' @export 
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' root <- root.node(g);
#' S.hier <- TPR.DAG(S, g, root, positive="children", bottomup="threshold.free", topdown="HTD", 
#' t=0, w=0, W=NULL, parallel=FALSE, ncores=1);
TPR.DAG <- function(S, g, root="00", positive="children", bottomup="threshold.free", topdown="HTD",
	t=0, w=0, W=NULL, parallel=FALSE, ncores=1){
	
	## Setting Check
	if(positive!="children" && positive!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && 
		bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau" || topdown!="HTD" && topdown!="GPAV")
		stop("TPR-DAG: positive or bottomup or topdown value misspelled", call.=FALSE);
	if(positive=="children" && bottomup=="tau")
		stop("TPR-DAG: tau is a descendants variants. Please set positive to descendants", call.=FALSE);
	if(bottomup=="threshold" || bottomup=="tau")
		w <- 0;
	if(bottomup=="threshold.free"){
		t <- 0;
		w <- 0;
	}
	if(bottomup=="weighted.threshold.free")
		t <-0;
	if(t==1 || w==1)
		warning("TPR-DAG: when t or w is equal to 1, TPR-DAG is reduced to HTD-DAG", call.=FALSE);	
	if(topdown=="GPAV" && parallel==TRUE && ncores<2)
		warning("GPAV: set ncores greater than 2 to exploit the GPAV parallel version", call.=FALSE);
	if(topdown=="GPAV" && parallel==FALSE && ncores>=2)
		warning("TPR-DAG: no GPAV parallel version is running, but ncores is higher or equal to 2.", 
			" Set 'ncores' to 1 to run the sequential version or set 'parallel' to TRUE to run the parallel version", call.=FALSE);
	if(topdown=="HTD" && (parallel==TRUE || ncores>=2))
		warning("TPR-DAG: does not exist a parallel version of HTD. The parameters 'parallel' and 'ncores' do not effect on 'HTD'", 
			". Set 'parallel' to FALSE and/or 'ncores' to 1 to avoid this warning message", call.=FALSE);	

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
		stop("TPR-DAG: mismatch between the number of nodes of the graph and the number of class of the scores matrix", call.=FALSE);

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
					children <- as.matrix(S[,chd.bup[[i]]]);	# genes considering children node 
					desc <-  as.matrix(S[,delta]);		# genes considering descendants nodes without children
				}else{
					desc <- as.matrix(S[,tmp]);
				}
				for(j in 1:length(parent)){
					if(bottomup=="threshold"){
						desc.set <- desc[j,] > t;    # positive descendants selection
						desc.pos <- desc[j,][desc.set];
						parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
					}else if(bottomup=="threshold.free"){
						desc.set <- desc[j,] > parent[j];	# positive descendants selection
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
						desc.set <- desc[j,] > parent[j];			# positive descendants (without children) selection
						desc.pos <- desc[j,][desc.set];
						child.set <- children[j,] > parent[j];  	# positive children selection
						child.pos <- children[j,][child.set];
						parent[j] <- t * ((parent[j] + sum(child.pos))/(1+length(child.pos))) + (1-t) * ((parent[j] + sum(desc.pos))/(1+length(desc.pos)));
					}
					S[,names(desc.bup[i])] <- parent;
				}
			}
		}
	}
	# top-down visit
	if(topdown=="HTD"){
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
	}else if(topdown=="GPAV"){
		if(parallel){
			S <- GPAV.parallel(S, g, W=W, ncores=ncores);
		}else{
			S <- GPAV.over.examples(S, W=W, g);
		}
	}	
	S <- S[,-which(colnames(S)==root)];
	return(S);
}


#' @name TPR-DAG-cross-validation
#' @title TPR-DAG cross-validation experiments
#' @seealso \code{\link{TPR-DAG-variants}}
#' @description High level function to correct the computed scores in a hierarchy according to the chosen ensemble algorithm 
#' @details The parametric hierarchical ensemble variants are cross-validated by maximizing in according to the metric
#' chosen in the parameter \code{metric}, that is F-measure (\code{\link{Multilabel.F.measure}}) or AUPRC (\code{\link{AUPRC}}).
#' @details The function checks if the number of classes between the flat scores matrix and the annotations matrix mismatched.
#' If so, the number of terms of the annotations matrix is shrunk to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well. N.B.: it is supposed that all the nodes of the subgraph are accessible from the root.
#' @param threshold range of threshold values to be tested in order to find the best threshold (\code{def:} \code{from:0.1}, 
#' \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but obviously the execution time will be higher.
#' Use this parameter only for the \emph{thresholded} variants; for the \emph{threshold-free} variant the functions automatically 
#' sets \code{threshold} to zero
#' @param weight range of weight values to be tested in order to find the best weight (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but obviously the execution time will be higher.
#' Use this parameter only for the \emph{weighted} variants; for the other variants the function automatically sets \code{weight} to zero
#' @param kk number of folds of the cross validation (\code{def: kk=5}) on which tuning the parameters \code{threshold}, \code{weight} and 
#' \code{tau} of the parametric variants of the hierarchical ensemble algorithms. For the non-parametric variants
#' (i.e. if \code{bottomup = threshold.free}) the function automatically sets \code{kk=NULL} if the input \code{kk}\eqn{\ne}\code{NULL} 
#' @param folds number of folds of the cross validation on which computing the performance metrics averaged across folds (\code{def. 5}).
#' If \code{folds=NULL}, the performance metrics are computed one-shot, otherwise the performance metrics are averaged across folds.
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. The parameter \code{seed} controls both the parameter \code{kk} and the parameter \code{folds}.
#' @param norm boolean value: should the flat scores matrix be normalized?
#' \itemize{
#' \item \code{TRUE} (\code{def.}): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} to set the on the fly 
#'	normalization method to apply among those possible.
#' }
#' @param norm.type can be one of the following three values: 
#'  \enumerate{
#'  \item \code{NULL} (\code{def.}): set \code{norm.type} to \code{NULL} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used. 
#'  }
#' @param positive choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#' 	\item \code{children} (\code{def.}): for each node are considered its positive children;
#' 	\item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#' 	\item \code{threshold.free} (\code{def.}): positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#' 	\item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#' 	\item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#' 	\item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#' 	\item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy. 
#'	NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positive=descendants};
#' }
#' @param topdown strategy to make the scores hierarchy-consistent. It can be one of the following values:
#' \itemize{
#' 	\item \code{HTD} (\code{def.}): \code{HTD-DAG} strategy is applied (\code{\link{HTD-DAG}});
#' 	\item \code{GPAV}: \code{GPAV} strategy is applied (\code{\link{GPAV}}).
#' }
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (by \code{def. W=NULL}), \code{W} is a unitary 
#' vector of the same length of the columns' number of the flat scores matrix (root node included). Set \code{W} only if \code{topdown=GPAV}.
#' @param parallel boolean value:
#' \itemize{
#'	\item \code{TRUE}: execute the parallel implementation of GPAV (\code{\link{GPAV.parallel}});
#'	\item \code{FALSE} (\code{def.}): execute the sequential implementation of GPAV (\code{\link{GPAV.over.examples}}).
#' }
#' Use \code{parallel} if and only if \code{topdown=GPAV}; otherwise set \code{parallel=FALSE}.
#' @param ncores number of cores to use for parallel execution (\code{def. 8}). Set \code{ncores=1} if \code{parallel=FALSE}, 
#' otherwise set \code{ncores} to the desired number of cores. Use \code{ncores} only if \code{topdown=GPAV}; otherwise set \code{parallel=1}.
#' @param recall.levels a vector with the desired recall levels (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}) to compute the 
#' the Precision at fixed Recall level (PXR)
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (\code{def.}): corresponds to the harmonic mean between the average precision and recall
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples
#' }
#' @param metric a string character specifying the performance metric on which to maximize the parametric ensemble variant. 
#' It can be one of the following values:
#' \enumerate{
#' \item \code{PRC}: the parametric ensemble variant is maximized on the basis of AUPRC (\code{\link{AUPRC}});
#' \item \code{FMAX}: the parametric ensemble variant is maximized on the basis of Fmax (\code{\link{Multilabel.F.measure}};
#' \item \code{NULL}: on the \code{threshold.free} variant none parameter optimization is needed, since the variant is non-parametric. 
#' So, if \code{bottomup=threshold.free} set \code{metric=NULL} (\code{def.}).
#' }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored
#' @param perf.dir relative path where the performance measures must be stored
#' @return Two \code{rda} files stored in the respective output directories:
#' \enumerate{
#' 	\item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' 	for each example and for each considered class. It is stored in the \code{hierScore.dir} directory.
#' 	\item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performance results:
#' 	\enumerate{
#' 		\item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'		\item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#' 		\item PXR results computed though \code{precision.at.given.recall.levels.over.classes} (\code{\link{PXR}});
#' 		\item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
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
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' threshold <- weight <- 0;
#' norm.type <- "MaxNorm";
#' positive <- "children";
#' bottomup <- "threshold.free";
#' topdown <- "HTD";
#' recall.levels <- seq(from=0.25, to=1, by=0.25);
#' Do.TPR.DAG(threshold=threshold, weight=weight, kk=NULL, folds=NULL, seed=NULL, norm=FALSE, 
#' norm.type=norm.type, positive=positive, bottomup=bottomup, topdown=topdown, W=NULL, 
#' parallel=FALSE, ncores=1, n.round=3, f.criterion="F", metric=NULL, recall.levels=recall.levels, 
#' flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
#' ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.TPR.DAG <- function(threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=0.9, by=0.1), 
	kk=5, folds=5, seed=23, norm=TRUE, norm.type=NULL, positive="children", bottomup="threshold.free", 
	topdown="HTD", W=NULL, parallel=FALSE, ncores=1, recall.levels=seq(from=0.1, to=1, by=0.1), n.round=3, 
	f.criterion="F", metric=NULL, flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
	flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	## Setting Check
	if(positive!="children" && positive!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && 
		bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau" || topdown!="HTD" && topdown!="GPAV")
		stop("TPR-DAG: positive or bottomup or topdown value misspelled", call.=FALSE);
	if(positive=="children" && bottomup=="tau")
		stop("TPR-DAG: tau is a descendants variants. Please set positive to descendants", call.=FALSE);
	if(bottomup=="threshold" || bottomup=="tau")
		weight <- 0;
	if(bottomup=="threshold.free"){
		threshold <- 0; 
		weight <- 0;
	}
	if(bottomup=="weighted.threshold.free")
		threshold <- 0;
	if(norm==FALSE && is.null(norm.type))
		stop("TPR-DAG: If norm is set to FALSE, you need to specify a normalization method among those available", call.=FALSE);
	if(norm==TRUE && !is.null(norm.type))
		warning("TPR-DAG: If norm is set to TRUE, the input flat matrix is already normalized.", 
			paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
	if((is.null(kk) || kk<=1) && bottomup!="threshold.free")
		stop("TPR-DAG: Smallest number of folds to define test and training set is 2. Set kk larger or equal to 2", call.=FALSE);
	if(!is.null(kk) && bottomup=="threshold.free")
		kk <- NULL;	
	if(f.criterion!="F" && f.criterion!="avF")
		stop("TPR-DAG: value of parameter 'f.criterion' misspelled");	
	if(metric!="FMAX" && metric!="PRC" && !is.null(metric))
		stop("TPR-DAG: value of parameter 'metric' misspelled", call.=FALSE);
	if(is.null(metric) && bottomup!="threshold.free")
		stop(paste0("TPR-DAG: metric cannot be NULL. The bottom-up approach '", bottomup, "' is parametric"), 
			". Please select the metric on which maximize according to those available", call.=FALSE); 
	if(!is.null(metric) && bottomup=="threshold.free") 
		warning("TPR-DAG: the bottom-up approach 'threshold.free' is non-parametric.", 
			paste0(" The chosen parameter metric '", metric, "' does not effect on 'threshold.free' bottom-up approach"), 
			". Set metric to NULL to avoid this warning message", call.=FALSE);
	if(is.null(seed) && bottomup!="threshold.free")
		warning("TPR-DAG: folds are generate without seed initialization", call.=FALSE);

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
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}else{
		S <- get(load(flat.path));
		S <- scores.normalization(norm.type=norm.type, S);
		cat(norm.type, "NORMALIZATION: DONE", "\n");
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}

	## check if |flat matrix classes| = |annotation matrix classes| 
	## if not the classes of annotation matrix are shrunk to those of flat matrix
	class.check <- ncol(S)!=ncol(ann);
	if(class.check){
		ann <- ann[,colnames(S)];
		nd <- c(root, colnames(S));
		g <- do.subgraph(nd, g, edgemode="directed");
	}

	## Compute FLAT PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
	PRC.flat <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
	AUC.flat <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
	PXR.flat <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
	FMM.flat <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
		b.per.example=TRUE, folds=folds, seed=seed);
	cat("FLAT PERFORMANCE: DONE", "\n");

	## Hierarchical Correction 
	if(bottomup=="threshold.free"){
		S <- TPR.DAG(S, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, t=0, w=0, W=W, 
			parallel=parallel, ncores=ncores);
		cat("HIERARCHICAL CORRECTION: DONE", "\n");
		## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated
		PRC.hier <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
		AUC.hier <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
		PXR.hier <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
		FMM.hier <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
			b.per.example=TRUE, folds=folds, seed=seed);
		cat("HIERARCHICAL PERFORMANCE: DONE", "\n");
		S.hier <- S;
		rm(S);
	}else{
		## Let's start k-fold crossing validation for choosing best threshold and weight maximizing on the selected metric
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
					pred.training <- TPR.DAG(training, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, 
						w=w, t=t, W=W, parallel=parallel, ncores=ncores);
					if(metric=="FMAX"){
						training.metric <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, 
							verbose=FALSE, b.per.example=FALSE)[["F"]];
					}else{
						training.metric <- AUPRC.single.over.classes(target.training, pred.training, folds=NULL, seed=NULL)$average;
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
							cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best threshold:",bestT, 
								"best weight:", bestW, sep="\t", "\n");
						}
					}
				}
			}
			## test set 
			pred.test <- TPR.DAG(test, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, t=bestT, w=bestW, 
				W=W, parallel=parallel, ncores=ncores);
			## assembling the hierarchical scores of each k sub-matrix
			S.hier <- rbind(S.hier, pred.test); 
		}
		## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the beginning matrix
		S.hier <- S.hier[rownames(S),];
		cat("HIERARCHICAL CORRECTION: DONE", "\n");
		## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated
		PRC.hier <- AUPRC.single.over.classes(ann, S.hier, folds=folds, seed=seed);
		AUC.hier <- AUROC.single.over.classes(ann, S.hier, folds=folds, seed=seed);
		PXR.hier <- precision.at.given.recall.levels.over.classes(ann, S, folds=folds, seed=seed, recall.levels=recall.levels);
		FMM.hier <- compute.Fmeasure.multilabel(ann, S.hier, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
			b.per.example=TRUE, folds=folds, seed=seed);
		cat("HIERARCHICAL PERFORMANCE: DONE", "\n");
		## remove no longer useful variables..
		rm(S, testIndex, pred.test, test, training, target.test, target.training); gc();
	}

	## Storing Results 
	if(positive=="children" && bottomup=="threshold.free" && topdown=="HTD")
		meth.name <- "tprTF";
	if(positive=="children" && bottomup=="threshold" && topdown=="HTD")
		meth.name <- "tprT";
	if(positive=="children" && bottomup=="weighted.threshold.free" && topdown=="HTD")
		meth.name <- "tprW";
	if(positive=="children" && bottomup=="weighted.threshold" && topdown=="HTD")
		meth.name <- "tprWT";
	if(positive=="descendants" && bottomup=="threshold.free" && topdown=="HTD")
		meth.name <- "descensTF";
	if(positive=="descendants" && bottomup=="threshold" && topdown=="HTD")
		meth.name <- "descensT";
	if(positive=="descendants" && bottomup=="weighted.threshold.free" && topdown=="HTD")
		meth.name <- "descensW";
	if(positive=="descendants" && bottomup=="weighted.threshold" && topdown=="HTD")
		meth.name <- "descensWT";
	if(positive=="descendants" && bottomup=="tau" && topdown=="HTD")
		meth.name <- "descensTAU";

	if(positive=="children" && bottomup=="threshold.free" && topdown=="GPAV")
		meth.name <- "ISOtprTF";
	if(positive=="children" && bottomup=="threshold" && topdown=="GPAV")
		meth.name <- "ISOtprT";
	if(positive=="children" && bottomup=="weighted.threshold.free" && topdown=="GPAV")
		meth.name <- "ISOtprW";
	if(positive=="children" && bottomup=="weighted.threshold" && topdown=="GPAV")
		meth.name <- "ISOtprWT";
	if(positive=="descendants" && bottomup=="threshold.free" && topdown=="GPAV")
		meth.name <- "ISOdescensTF";
	if(positive=="descendants" && bottomup=="threshold" && topdown=="GPAV")
		meth.name <- "ISOdescensT";
	if(positive=="descendants" && bottomup=="weighted.threshold.free" && topdown=="GPAV")
		meth.name <- "ISOdescensW";
	if(positive=="descendants" && bottomup=="weighted.threshold" && topdown=="GPAV")
		meth.name <- "ISOdescensWT";
	if(positive=="descendants" && bottomup=="tau" && topdown=="GPAV")
		meth.name <- "ISOdescensTAU";

	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier, 
			file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier,
			file=paste0(perf.dir, "PerfMeas.", norm.type, ".", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
	}
}

#' @name TPR-DAG-holdout
#' @title TPR-DAG holdout experiments
#' @seealso \code{\link{TPR-DAG-variants}}
#' @description High level function to correct the computed scores in a hierarchy according to the chosen ensemble algorithm 
#' @details The parametric hierarchical ensemble variants are cross-validated by maximizing in according to the metric
#' chosen in the parameter \code{metric}, that is F-measure (\code{\link{Multilabel.F.measure}}) or AUPRC (\code{\link{AUPRC}}).
#' @details The function checks if the number of classes between the flat scores matrix and the annotations matrix mismatched.
#' If so, the number of terms of the annotations matrix is shrunk to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well. N.B.: it is supposed that all the nodes of the subgraph are accessible from the root.
#' @param threshold range of threshold values to be tested in order to find the best threshold (\code{def:} \code{from:0.1}, 
#' \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but obviously the execution time will be higher.
#' Use this parameter only for the \emph{thresholded} variants; for the \emph{threshold-free} variant the functions automatically 
#' sets \code{threshold} to zero
#' @param weight range of weight values to be tested in order to find the best weight (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best threshold is, but obviously the execution time will be higher.
#' Use this parameter only for the \emph{weighted} variants; for the other variants the function automatically sets \code{weight} to zero
#' @param kk number of folds of the cross validation (\code{def: kk=5}) on which tuning the parameters \code{threshold}, \code{weight} and 
#' \code{tau} of the parametric variants of the hierarchical ensemble algorithms. For the non-parametric variants
#' (i.e. if \code{bottomup = threshold.free}) the function automatically sets \code{kk=NULL} if the input \code{kk}\eqn{\ne}\code{NULL} 
#' @param folds number of folds of the cross validation on which computing the performance metrics averaged across folds (\code{def. 5}).
#' If \code{folds=NULL}, the performance metrics are computed one-shot, otherwise the performance metrics are averaged across folds.
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. The parameter \code{seed} controls both the parameter \code{kk} and the parameter \code{folds}.
#' @param norm boolean value: should the flat scores matrix be normalized?
#' \itemize{
#' \item \code{TRUE} (\code{def.}): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} to set the on the fly 
#'	normalization method to apply among those possible.
#' }
#' @param norm.type can be one of the following three values: 
#'  \enumerate{
#'  \item \code{NULL} (\code{def.}): set \code{norm.type} to \code{NULL} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used. 
#'  }
#' @param positive choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#' 	\item \code{children} (\code{def.}): for each node are considered its positive children;
#' 	\item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#' 	\item \code{threshold.free} (\code{def.}): positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#' 	\item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#' 	\item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#' 	\item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#' 	\item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy. 
#'	NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positive=descendants};
#' }
#' @param topdown strategy to make the scores hierarchy-consistent. It can be one of the following values:
#' \itemize{
#' 	\item \code{HTD} (\code{def.}): \code{HTD-DAG} strategy is applied (\code{\link{HTD-DAG}});
#' 	\item \code{GPAV}: \code{GPAV} strategy is applied (\code{\link{GPAV}}).
#' }
#' @param W vector of weight relative to a single example. If the vector \code{W} is not specified (by \code{def.} \code{W=NULL}), \code{W} is a unitary 
#' vector of the same length of the columns' number of the flat scores matrix (root node included). Set \code{W} only if \code{topdown=GPAV}.
#' @param parallel boolean value:
#' \itemize{
#'	\item \code{TRUE}: execute the parallel implementation of GPAV (\code{\link{GPAV.parallel}});
#'	\item \code{FALSE} (\code{def.}): execute the sequential implementation of GPAV (\code{\link{GPAV.over.examples}}).
#' }
#' Use \code{parallel} if and only if \code{topdown=GPAV}; otherwise set \code{parallel=FALSE}.
#' @param ncores number of cores to use for parallel execution (\code{def. 8}). Set \code{ncores=1} if \code{parallel=FALSE}, 
#' otherwise set \code{ncores} to the desired number of cores.
#' Use \code{ncores} if and only if \code{topdown=GPAV}; otherwise set \code{parallel=1}.
#' @param recall.levels a vector with the desired recall levels (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}) to compute the 
#' the Precision at fixed Recall level (PXR)
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (def.): corresponds to the harmonic mean between the average precision and recall
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples
#' }
#' @param metric a string character specifying the performance metric on which to maximize the parametric ensemble variant. 
#' It can be one of the following values:
#' \enumerate{
#' \item \code{PRC}: the parametric ensemble variant is maximized on the basis of AUPRC (\code{\link{AUPRC}});
#' \item \code{FMAX}: the parametric ensemble variant is maximized on the basis of Fmax (\code{\link{Multilabel.F.measure}});
#' \item \code{NULL}: on the \code{threshold.free} variant none parameter optimization is needed, since the variant is non-parametric. 
#' So, if \code{bottomup=threshold.free} set \code{metric=NULL} (\code{def.}).
#' }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param ind.test.set name of the file containing a vector of integer numbers corresponding to the indices of the elements (rows) of scores 
#' matrix to be used in the	test set 
#' @param ind.dir relative path to folder where \code{ind.test.set} is stored
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored
#' @param perf.dir relative path where the performance measures must be stored
#' @return Two \code{rda} files stored in the respective output directories:
#' \enumerate{
#' 	\item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' 	for each example and for each considered class. It is stored in the \code{hierScore.dir} directory.
#' 	\item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performace results:
#' 	\enumerate{
#' 		\item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'		\item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#' 		\item PXR results computed though \code{precision.at.given.recall.levels.over.classes} (\code{\link{PXR}});
#' 		\item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
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
#' ind.test.set <- "test.index";
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' threshold <- weight <- 0;
#' norm.type <- "MaxNorm";
#' positive <- "children";
#' bottomup <- "threshold.free";
#' topdown <- "HTD";
#' recall.levels <- seq(from=0.25, to=1, by=0.25);
#' Do.TPR.DAG.holdout(threshold=threshold, weight=weight, kk=NULL, folds=NULL, seed=NULL, norm=FALSE, 
#' norm.type=norm.type, positive=positive, bottomup=bottomup, topdown=topdown, W=NULL, 
#' parallel=FALSE, ncores=1, recall.levels=recall.levels, n.round=3, f.criterion="F", metric=NULL,
#' flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, 
#' ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.TPR.DAG.holdout <- function(threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=1, by=0.1), kk=5, 
	folds=5, seed=23, norm=TRUE, norm.type=NULL, positive="children", bottomup="threshold.free", topdown="HTD",
	W=NULL, parallel=FALSE, ncores=1, recall.levels=seq(from=0.1, to=1, by=0.1), n.round=3, f.criterion="F", metric=NULL,
	flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, 
	flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir){

	## Setting Check
	if(positive!="children" && positive!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && 
		bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau" || topdown!="HTD" && topdown!="GPAV")
		stop("TPR-DAG: positive or bottomup or topdown value misspelled", call.=FALSE);
	if(positive=="children" && bottomup=="tau")
		stop("TPR-DAG: tau is a descendants variants. Please set positive to descendants", call.=FALSE);
	if(bottomup=="threshold" || bottomup=="tau")weight<-0;
	if(bottomup=="threshold.free"){
		threshold <- 0; 
		weight <- 0;
	}
	if(bottomup=="weighted.threshold.free")
		threshold<-0;
	if(norm==FALSE && is.null(norm.type))
		stop("TPR-DAG: If norm is set to FALSE, you need to specify a normalization method among those available", call.=FALSE);
	if(f.criterion!="F" && f.criterion!="avF")
		stop("TPR-DAG: value of parameter 'f.criterion' misspelled", call.=FALSE);	
	if(metric!="FMAX" && metric!="PRC" && !is.null(metric))
		stop("TPR-DAG: value of parameter 'metric' misspelled", call.=FALSE);
	if(norm==TRUE && !is.null(norm.type))
		warning("TPR-DAG: If norm is set to TRUE, the input flat matrix is already normalized.", 
			paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
	if((is.null(kk) || kk<=1) && bottomup!="threshold.free")
		stop("TPR-DAG: Smallest number of folds to define test and training set is 2. Set kk larger or equal to 2", call.=FALSE);
	if(!is.null(kk) && bottomup=="threshold.free")
		kk <- NULL;			
	if(is.null(metric) && bottomup!="threshold.free")
		stop(paste0("TPR-DAG: metric cannot be NULL. The bottom-up approach '", bottomup, "' is parametric"), 
			". Please select the metric on which maximize according to those available", call.=FALSE); 
	if(!is.null(metric) && bottomup=="threshold.free") 
		warning("TPR-DAG: the bottom-up approach 'threshold.free' is non-parametric.", 
			paste0(" The chosen parameter metric '", metric, "' does not effect on 'threshold.free' bottom-up approach"), 
			". Set metric to NULL to avoid this warning message", call.=FALSE);
	if(is.null(seed) && bottomup!="threshold.free")
		warning("TPR-DAG: folds are generate without seed initialization", call.=FALSE);

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
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}else{
		S <- get(load(flat.path));
		S <- scores.normalization(norm.type=norm.type, S);
		cat(norm.type, "NORMALIZATION: DONE", "\n");
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}

	## check if |flat matrix classes| = |annotation matrix classes| 
	## if not the classes of annotation matrix are shrunk to those of flat matrix
	class.check <- ncol(S)!=ncol(ann);
	if(class.check){
		ann <- ann[,colnames(S)];
		nd <- c(root, colnames(S));
		g <- do.subgraph(nd, g, edgemode="directed");
	}

	## scores flat matrix are shrunk to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];

	## annotation table are shrunk to test and training test respectively.
	ann.test <- ann[ind.test,];
	ann.training <- ann[-ind.test,];

	## Compute FLAT PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test, folds=folds, seed=seed);
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test, folds=folds, seed=seed);
	PXR.flat <- precision.at.given.recall.levels.over.classes(ann.test, S.test, folds=folds, seed=seed, recall.levels=recall.levels);
	FMM.flat <- compute.Fmeasure.multilabel(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE,
		b.per.example=TRUE, folds=folds, seed=seed);
	cat("FLAT PERFORMANCE: DONE", "\n");

	## Hierarchical Correction 
	if(bottomup=="threshold.free"){
		S.test <- TPR.DAG(S.test, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, 
			t=0, w=0, W=W, parallel=parallel, ncores=ncores);
		cat("HIERARCHICAL CORRECTION: DONE", "\n");
		## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
		PRC.hier <- AUPRC.single.over.classes(ann.test, S.test, folds=folds, seed=seed);
		AUC.hier <- AUROC.single.over.classes(ann.test, S.test, folds=folds, seed=seed);
		PXR.hier <- precision.at.given.recall.levels.over.classes(ann.test, S.test, folds=folds, seed=seed, recall.levels=recall.levels);
		FMM.hier <- compute.Fmeasure.multilabel(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
			b.per.example=TRUE, folds=folds, seed=seed);
		cat("HIERARCHICAL PERFORMANCE: DONE", "\n");
		S.hier <- S.test;
		rm(S, S.test);
	}else{
		## Let's start k-fold crossing validation for choosing best threshold and weight maximizing on the selected metric
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
					pred.training <- TPR.DAG(training, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, 
						w=w, t=t, W=W, parallel=parallel, ncores=ncores);
					if(metric=="FMAX"){
						training.metric <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, 
							verbose=FALSE, b.per.example=FALSE)[["F"]];
					}else{
						training.metric <- AUPRC.single.over.classes(target.training, pred.training, folds=NULL, seed=NULL)$average;	 
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
							cat("training fold:", k, paste0("top ", metric," avg found:"), top.metric, "best threshold:",bestT, 
								"best weight:", bestW, sep="\t", "\n");
						}
					}
				}
			}
		}
		S.test <- TPR.DAG(S.test, g, root=root, positive=positive, bottomup=bottomup, topdown=topdown, t=bestT, w=bestW, 
			W=W, parallel=parallel, ncores=ncores);
		cat("HIERARCHICAL CORRECTION: DONE", "\n");

		## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
		PRC.hier <- AUPRC.single.over.classes(ann.test, S.test, folds=folds, seed=seed);
		AUC.hier <- AUROC.single.over.classes(ann.test, S.test, folds=folds, seed=seed);
		PXR.hier <- precision.at.given.recall.levels.over.classes(ann.test, S.test, folds=folds, seed=seed, recall.levels=recall.levels);
		FMM.hier <- compute.Fmeasure.multilabel(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
			b.per.example=TRUE, folds=folds, seed=seed);
		cat("HIERARCHICAL PERFORMANCE: DONE", "\n");
		
		## storing the hierarchical matrix
		S.hier <- S.test;
		rm(S, S.test, S.training, training);
	}
	## Storing Results 
	if(positive=="children" && bottomup=="threshold.free" && topdown=="HTD")
		meth.name <- "tprTF";
	if(positive=="children" && bottomup=="threshold" && topdown=="HTD")
		meth.name <- "tprT";
	if(positive=="children" && bottomup=="weighted.threshold.free" && topdown=="HTD")
		meth.name <- "tprW";
	if(positive=="children" && bottomup=="weighted.threshold" && topdown=="HTD")
		meth.name <- "tprWT";

	if(positive=="descendants" && bottomup=="threshold.free" && topdown=="HTD")
		meth.name <- "descensTF";
	if(positive=="descendants" && bottomup=="threshold" && topdown=="HTD")
		meth.name <- "descensT";
	if(positive=="descendants" && bottomup=="weighted.threshold.free" && topdown=="HTD")
		meth.name <- "descensW";
	if(positive=="descendants" && bottomup=="weighted.threshold" && topdown=="HTD")
		meth.name <- "descensWT";
	if(positive=="descendants" && bottomup=="tau" && topdown=="HTD")
		meth.name <- "descensTAU";

	if(positive=="children" && bottomup=="threshold.free" && topdown=="GPAV")
		meth.name <- "ISOtprTF";
	if(positive=="children" && bottomup=="threshold" && topdown=="GPAV")
		meth.name <- "ISOtprT";
	if(positive=="children" && bottomup=="weighted.threshold.free" && topdown=="GPAV")
		meth.name <- "ISOtprW";
	if(positive=="children" && bottomup=="weighted.threshold" && topdown=="GPAV")
		meth.name <- "ISOtprWT";
	
	if(positive=="descendants" && bottomup=="threshold.free" && topdown=="GPAV")
		meth.name <- "ISOdescensTF";
	if(positive=="descendants" && bottomup=="threshold" && topdown=="GPAV")
		meth.name <- "ISOdescensT";
	if(positive=="descendants" && bottomup=="weighted.threshold.free" && topdown=="GPAV")
		meth.name <- "ISOdescensW";
	if(positive=="descendants" && bottomup=="weighted.threshold" && topdown=="GPAV")
		meth.name <- "ISOdescensWT";
	if(positive=="descendants" && bottomup=="tau" && topdown=="GPAV")
		meth.name <- "ISOdescensTAU";

	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier,
			file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier,
			file=paste0(perf.dir, "PerfMeas.", norm.type, ".", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
	}
}
