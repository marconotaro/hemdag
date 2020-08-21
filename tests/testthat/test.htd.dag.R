library(HEMDAG);

context("test htd-dag");

check.graph <- function(){ 
    if(requireNamespace("graph", quietly=TRUE)){ 
        TRUE 
    }else{ 
        FALSE 
    } 
}

make.graph <- function(){
    if (!check.graph()){ skip("graph package cannot be loaded"); }

    V <- LETTERS[1:10];
    edL <- vector("list", length=length(V));
    names(edL) <- V;
    edL[[1]] <- list(edges=c(2,3,7,5));
    edL[[2]] <- list(edges=c(4,6));
    edL[[3]] <- list(edges=c(6,8)); 
    edL[[4]] <- list(edges=c(6,9));
    edL[[5]] <- list(edges=c(8));
    edL[[6]] <- list(edges=c(8,9,10));
    edL[[7]] <- list(edges=c(8,5));
    edL[[8]] <- list(edges=c());
    edL[[9]] <- list(edges=c());
    g <- graph::graphNEL(nodes=V, edgeL=edL, edgemode="directed");   
    return(g);
    ## optional: plotting the DAG g
    ## library(Rgraphviz); plot(g);
}

make.scores <- function(){
    pr1 <- c(0.9, 0.4, 0.4, 0.3, 0.6, 0.4, 0.5, 0.9, 0.3, 0.8);
    pr2 <- c(0.7, 0.9, 0.8, 0.9, 0.5, 0.2, 0.3, 0.5, 0.7, 1.0);
    pr3 <- c(0.1, 0.9, 0.7, 1.0, 0.9, 0.4, 0.8, 0.9, 0.1, 0.8);
    S <- rbind(pr1,pr2,pr3);
    colnames(S) <- LETTERS[1:length(pr1)];
    return(S);
}

make.ann <- function(){
    pr1 <- c(1,1,1,1,0,1,0,0,0,1);
    pr2 <- c(1,1,1,1,1,1,1,1,0,0);
    pr3 <- c(1,1,1,1,0,1,0,0,0,1);
    ann <- rbind(pr1,pr2,pr3);
    colnames(ann) <- c("A","B","C","D","E","F","G","H","I","J");
    return(ann);
}

test_that("htd works", {
    g <- make.graph();
    root <- root.node(g);
    
    ## test with root
    S  <- make.scores();
    S.htd  <- htd(S, g, root);
    pr1 <- c(0.9, 0.4, 0.4, 0.3, 0.5, 0.3, 0.5, 0.3, 0.3, 0.3);
    pr2 <- c(0.7, 0.7, 0.7, 0.7, 0.3, 0.2, 0.3, 0.2, 0.2, 0.2);
    pr3 <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
    S.check <- rbind(pr1,pr2,pr3);
    colnames(S.check) <- LETTERS[1:length(pr1)];
    expect_equal(S.htd, S.check);

    ## test without root
    S2 <- S[,-which(colnames(S) %in% root)]; 
    S.htd2 <- htd(S2, g, root);
    pr1 <- c(1, 0.4, 0.4, 0.3, 0.5, 0.3, 0.5, 0.3, 0.3, 0.3);
    pr2 <- c(1, 0.9, 0.8, 0.9, 0.3, 0.2, 0.3, 0.2, 0.2, 0.2);
    pr3 <- c(1, 0.9, 0.7, 0.9, 0.8, 0.4, 0.8, 0.4, 0.1, 0.4);
    S.check2 <- rbind(pr1,pr2,pr3);
    colnames(S.check2) <- LETTERS[1:length(pr1)];
    expect_equal(S.htd2, S.check2);

    ## test class mismatch
    S.error <- S[,-which(colnames(S) %in% c("D","H"))];
    expect_error(htd(S.error, g, root),
        "htd: the number of nodes of the graph and the number of classes of the scores matrix does not match");
})










