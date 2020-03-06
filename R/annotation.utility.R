################################################################
## Utility functions to process and analyze annotation tables ##
################################################################

#' @title specific annotation matrix
#' @description Construct the labels matrix of the most specific OBO terms (as GO or HPO).
#' @details The input plain text file representing the most specific associations gene-OBO term can be obtained by cloning the GitHub repository
#' \href{https://github.com/marconotaro/obogaf-parser}{obogaf-parser}, a perl5 module specifically designed to handle HPO and GO obo file and 
#' their gene annotation file (gaf file).
#' @param file text file representing the most specific associations gene-OBO. The file must be written as sequence of rows. 
#' Each row represents a gene/protein and all its associations with an ontology term pipe separated, \emph{e.g.: gene1 |obo1|...|oboN}.  
#' The input example file used here (\code{def: "gene2pheno.txt"}) shows the gene and all its associations with an HPO terms.
#' @return the annotation matrix of the most specific annotations (0/1): rows are genes and columns are HPO terms.
#' Let's denote \eqn{M} the labels matrix. If \eqn{M[i,j]=1}, means that the gene \eqn{i} is annotated with the class \eqn{j}, otherwise \eqn{M[i,j]=0}.
#' @export
#' @examples
#' gene2pheno <- system.file("extdata/gene2pheno.txt.gz", package="HEMDAG");
#' spec.ann <- specific.annotation.matrix(file=gene2pheno);
specific.annotation.matrix <- function(file="gene2pheno.txt.gz"){
    tmp <- strsplit(file, "[.,/,_]")[[1]];
    if(any(tmp %in% "gz")){
        con <- gzfile(file);
        line <- readLines(con);
        close(con);
    }else{
        line <- readLines(file);
    }
    tmp <- strsplit(line, split="[ |]"); #nb: space before pipe useful to separate gene from obo terms
    genenames <- c();
    for(i in 1:length(tmp)) genenames <- c(genenames,tmp[[i]][1]);
    ann.list <- list();
    for(i in 1:length(tmp)) ann.list[[i]] <- unique(tmp[[i]])[-1];
    names(ann.list) <- genenames;
    oboID <- unique(unlist(ann.list));
    n.genes <- length(genenames);
    n.oboID <- length(oboID);
    m <- matrix(integer(n.genes * n.oboID), nrow=n.genes);
    rownames(m) <- genenames;
    colnames(m) <- oboID;
    for (i in genenames){
        spec.ann <- ann.list[[i]]; 
        m[i, spec.ann] <- 1;  
    }
    charcheck <- any(suppressWarnings(is.na(as.numeric(genenames))));
    if(charcheck){
        m <- m[sort(rownames(m)),sort(colnames(m))];
    }else{
        rname <- as.character(sort(as.numeric(genenames)));
        m <- m[rname, sort(colnames(m))];
    }
    return(m);
}

#' @title Specific annotations list
#' @description Construct a list of the most specific annotations starting from the table of the most specific annotations.
#' @param ann annotation matrix (0/1). Rows are examples and columns are most specific terms. It must be a named matrix. 
#' @return a named list, where the names of each component correspond to an examples (genes) and the elements of each component 
#' are the most specific classes associated to that genes.
#' @seealso \code{\link{specific.annotation.matrix}}
#' @export
#' @examples
#' data(labels);
#' spec.list <- specific.annotation.list(L);
specific.annotation.list <- function(ann){
 ann.list <- apply(ann, 1, function(gene){
        terms <- which(gene==1);
        return(names(gene[terms]));
    });
    return(ann.list);
}

#' @title Transitive closure of annotations 
#' @description Performs the transitive closure of the annotations using ancestors and the most specific annotation table.
#' The annotations are propagated from bottom to top, enriching the most specific annotations table.
#' The rows of the matrix correspond to the genes of the most specific annotation table and the columns to the OBO terms/classes.
#' @param ann.spec the annotation matrix of the most specific annotations (0/1): rows are genes and columns are OBO terms.
#' @param anc list of the ancestors of the ontology. 
#' @return an annotation table T: rows correspond to genes and columns to OBO terms. \eqn{T[i,j]=1} means that gene \eqn{i} is annotated for the term \eqn{j},
#' \eqn{T[i,j]=0} means that gene \eqn{i} is not annotated for the term \eqn{j}.
#' @seealso \code{\link{specific.annotation.matrix}}, \code{\link{build.ancestors}}
#' @export
#' @examples
#' data(graph);
#' data(labels);
#' anc <- build.ancestors(g);
#' tca <- transitive.closure.annotations(L, anc);
transitive.closure.annotations <- function(ann.spec, anc){
    ## costructiion of annotation list
    ann.list <- specific.annotation.list(ann.spec);
    ## cotruction the full empty annotation matrix
    genes <- rownames(ann.spec);
    n.genes <- length(genes);
    oboIDs <- names(anc);
    n.oboID <- length(anc);
    obo.ann <- matrix(numeric(n.oboID * n.genes), nrow=n.genes, ncol=n.oboID);    #empty label matrix
    dimnames(obo.ann) <- list(genes,oboIDs);
    ## fill the full empty annotation matrix with the most specific annotation     
    obo.spec.term <- colnames(ann.spec); # the most specific OBO terms
    # might happen that there are same OBO IDs that are classified as "obsolete" in obo file, but that still exist in the annotation file 
    obo.spec.term.sel <- oboIDs[oboIDs  %in% obo.spec.term]; # removing obsolete OBO terms...
    obo.ann[genes,obo.spec.term.sel] <- ann.spec[,obo.spec.term.sel];    
    ## transitive closure: annotation propagation from the most specific nodes to all its ancestors
    for (i in genes){
        spec.ann <- ann.list[[i]];
        all.anc <- lapply(spec.ann, function(x) return(anc[[x]]));
        all.anc <- unique(unlist(all.anc));
        obo.ann[i, all.anc] <- 1;  # setting the annotations derived by transitive closure
    }
    ## remove OBO empty terms 
    obo.ann <- obo.ann[,colSums(obo.ann)!=0];
    return(obo.ann);
}

#' @title Full annotation matrix
#' @description Construct a full annotations table using ancestors and the most specific annotations table w.r.t. a given weighted adjacency matrix (wadj). 
#' The rows of the full annotation matrix correspond to all the examples of the given weighted adjacency matrix and the columns to the class/terms.
#' The transitive closure of the annotations is performed. 
#' @details The examples present in the annotation matrix (\code{ann.spec}) but not in the adjacency weighted matrix (\code{W}) are purged.
#' @param W symmetric adjacency weighted matrix of the graph. 
#' @param anc list of the ancestors of the ontology. 
#' @param ann.spec the annotation matrix of the most specific annotations (0/1): rows are genes and columns are classes.
#' @return a full annotation table T, that is a matrix in which the transitive closure of annotations was performed. 
#' Rows correspond to genes of the weighted adjacency matrix and columns to terms. 
#' \eqn{T[i,j]=1} means that gene \eqn{i} is annotated for the term \eqn{j}, \eqn{T[i,j]=0} means that gene \eqn{i} is not annotated for the term \eqn{j}.
#' @seealso \code{\link{weighted.adjacency.matrix}}, \code{\link{build.ancestors}}, \cr
#' \code{\link{specific.annotation.matrix}}, \code{\link{transitive.closure.annotations}}
#' @export
#' @examples
#' data(wadj);
#' data(graph);
#' data(labels);
#' anc <- build.ancestors(g);
#' full.ann <- full.annotation.matrix(W, anc, L);
full.annotation.matrix <- function(W, anc, ann.spec){
    ## construction of annotation list
    ann.list <- specific.annotation.list(ann.spec);
    ## construction the full empty annotation matrix
    genes <- rownames(W);
    n.genes <- length(genes);
    oboIDs <- names(anc);
    n.oboID <- length(anc);
    obo.ann <- matrix(numeric(n.oboID * n.genes), nrow=n.genes, ncol=n.oboID);    #empty label matrix
    dimnames(obo.ann) <- list(genes,oboIDs);
    ## fill the full empty annotation matrix with the most specific annotation 
    genes2obo <- rownames(ann.spec);              # all genes that are associated with OBO terms
    genes.sel <- genes[genes %in% genes2obo];     # genes 2 OBO terms 2 entrez id of wadj
    obo.spec.term <- colnames(ann.spec);          # the most specific OBO terms
    #might happen that there are same obo IDs that are classified as "obsolete" in obo file, but that still exist in the annotation file (e.g. build 1233)
    obo.spec.term.sel <- oboIDs[oboIDs  %in% obo.spec.term]; # removing obsolete OBO terms...
    obo.ann[genes.sel,obo.spec.term.sel] <- ann.spec[genes.sel,obo.spec.term.sel];    # setting the most specific annotations
    ## transitive closure: annotation propagation from the most specific nodes to all its ancestors
    for (i in genes){
        spec.ann <- ann.list[[i]];
        all.anc <- lapply(spec.ann, function(x) return(anc[[x]]));
        all.anc <- unique(unlist(all.anc));
        obo.ann[i, all.anc] <- 1;  # setting the annotations derived by transitive closure
    }
    ## remove OBO empty terms 
    obo.ann <- obo.ann[,colSums(obo.ann)!=0];
    return(obo.ann);
}

#' @title Do full annotation matrix
#' @description High-level function to obtain a full annotation matrix, that is a matrix in which the transitive closure of annotations was performed, 
#' respect to a given weighted adjacency matrix.
#' @param anc.file.name name of the file containing the list for each node the list of all its ancestor (without \code{rda} extension).
#' @param anc.dir relative path to directory where the ancestor file is stored. 
#' @param net.file name of the file containing the weighted adjacency matrix of the graph (without \code{rda} extension).
#' @param net.dir relative path to directory where the weighted adjacency matrix is stored.
#' @param ann.file.name name of the file containing the matrix of the most specific annotations (without \code{rda} extension).
#' @param ann.dir relative path to directory where the matrix of the most specific annotation is stored.
#' @param output.name name of the output file without \code{rda} extension.
#' @param output.dir relative path to directory where the output file must be stored.   
#' @return a full annotation matrix T, that is a matrix in which the transitive closure of annotations was performed.
#' Rows correspond to genes of the input weighted adjacency matrix and columns to terms. 
#' \eqn{T[i,j]=1} means that gene \eqn{i} is annotated for the term \eqn{j}, \eqn{T[i,j]=0} means that gene \eqn{i} is not annotated for the term \eqn{j}.
#' @seealso \code{\link{full.annotation.matrix}}
#' @export
#' @examples
#' data(graph);
#' data(labels);
#' data(wadj);
#' anc <- build.ancestors(g);
#' tmpdir <- paste0(tempdir(),"/");
#' save(g, file=paste0(tmpdir,"graph.rda"));
#' save(L, file=paste0(tmpdir,"labels.rda"));
#' save(W, file=paste0(tmpdir,"wadj.rda"));
#' save(anc, file=paste0(tmpdir,"ancestors.rda"));
#' anc.dir <- net.dir <- ann.dir <- output.dir <- tmpdir;
#' anc.file.name <- "ancestors";
#' net.file <- "wadj";
#' ann.file.name <- "labels";
#' output.name <- "full.ann.matrix";
#' Do.full.annotation.matrix(anc.file.name=anc.file.name, anc.dir=anc.dir, net.file=net.file, 
#'    net.dir=net.dir, ann.file.name=ann.file.name, ann.dir=ann.dir, output.name=output.name, 
#'     output.dir=output.dir);
Do.full.annotation.matrix <- function(anc.file.name=anc.file.name, anc.dir=anc.dir, net.file=net.file, net.dir=net.dir, 
    ann.file.name=ann.file.name, ann.dir=ann.dir, output.name=output.name, output.dir=output.dir){
    ## loading list of ancestors
    anc.path <- paste0(anc.dir, anc.file.name, ".rda");
    anc.name <- load(anc.path);
    anc <- eval(parse(text=anc.name));  
    ## loading wadj matrix
    net.path <- paste0(net.dir, net.file, ".rda");
    net.name <- load(net.path);
    W <- eval(parse(text=net.name));  
    ## loading the specific annotation matrix
    ann.path <- paste0(ann.dir, ann.file.name, ".rda");
    ann.name <- load(ann.path);
    ann.spec <- eval(parse(text=ann.name));
    ## costruction of full OBO annotation matrix
    ann <- full.annotation.matrix(W=W, anc=anc, ann.spec=ann.spec);
    ## saving labels matrix
    ann.file <- paste0(output.dir, output.name, ".rda");
    save(ann, file=ann.file, compress=TRUE);
}

#' @title Build submatrix
#' @title Build an annotation matrix with only those terms having more than n annotations.
#' @description Terms having less than n annotations are pruned. Terms having exactly n annotations are discarded as well.
#' @param ann the annotation matrix (0/1). Rows are examples and columns are classes. 
#' @param n integer number of annotations to be pruned.
#' @return Matrix of annotations having only those terms with more than n annotations.
#' @export
#' @examples
#' data(labels);
#' subm <- do.submatrix(L,5);
do.submatrix <- function(ann,n){
    ann.sel <- ann[,colSums(ann)>n];
    return(ann.sel);
}

#' @title Annotation matrix checker
#' @description This function assess the integrity of an annotation table in which a transitive closure of annotations was performed.
#' @param anc list of the ancestors of the ontology. 
#' @param ann.spec the annotation matrix of the most specific annotations (0/1): rows are genes and columns are terms. 
#' @param ann the full annotation matrix (0/1), that is the matrix in which the transitive closure of the annotation was performed.
#' Rows are examples and columns are classes. 
#' @return If the transitive closure of the annotations is well performed "OK" is returned, otherwise a message error is printed on the stdout.
#' @seealso \code{\link{build.ancestors}}, \code{\link{transitive.closure.annotations}}, \code{\link{full.annotation.matrix}}
#' @export
#' @examples
#' data(graph);
#' data(labels);
#' anc <- build.ancestors(g);
#' tca <- transitive.closure.annotations(L, anc);
#' check.annotation.matrix.integrity(anc, L, tca);
check.annotation.matrix.integrity <- function(anc, ann.spec, ann){
    ## construction of annotation list
    ann.list <- specific.annotation.list(ann.spec);
    genes <- rownames(ann);
    check <- c();
    for (i in genes){
        spec.ann <- which(ann[i,]==1);
        len.ann <- length(spec.ann);
        all.anc <- lapply(ann.list[[i]], function(x) return(anc[[x]]));
        all.anc <- unique(unlist(all.anc));
        len.anc <- length(all.anc);
        cmp <- len.anc == len.ann;
        if(cmp==TRUE){
            check <- c(check,"OK");
        } else {
            check <- c(check,"NOTOK");
        }
    }
    names(check) <- genes;
    violated <- any(check!="OK");
    if(violated){
        n <- names(check)[check=="NOTOK"];
        cat("check.annotation.matrix: NOT_OK. Transitive closure NOT RESPECTED", "\n");
    }else{
        cat("check.annotation.matrix: OK", "\n");    
    }
}
