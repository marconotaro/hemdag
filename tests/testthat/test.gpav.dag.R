library(HEMDAG);

context("test gpav-dag");

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

test_that("adj.upper.tri works",{
    g <- make.graph();
    adj <- adj.upper.tri(g);

    H <- c(0,0,0,1,0,0,1,1,1,0)
    I <- c(0,0,0,1,1,0,0,0,0,0)
    J <- c(0,0,0,1,0,0,0,0,0,0)
    F <- c(0,0,0,0,1,1,1,0,0,0)
    D <- c(0,0,0,0,0,1,0,0,0,0)
    B <- c(0,0,0,0,0,0,0,0,0,1)
    C <- c(0,0,0,0,0,0,0,0,0,1)
    E <- c(0,0,0,0,0,0,0,0,1,1)
    G <- c(0,0,0,0,0,0,0,0,0,1)
    A <- c(0,0,0,0,0,0,0,0,0,0)
    adj.control <- rbind(H,I,J,F,D,B,C,E,G,A)
    colnames(adj.control) <- c("H","I","J","F","D","B","C","E","G","A");

    expect_equal(adj, adj.control);
})

test_that("gpav works", {
    g <- make.graph();
    adj <- adj.upper.tri(g);

    Y <- c(0.9, 0.4, 0.4, 0.3, 0.6, 0.4, 0.5, 0.9, 0.3, 0.8);
    names(Y) <- c("A","B","C","D","E","F","G","H","I","J");
    Y.trim <- c(0.9, 0.4, 0.4, 0.3, 0.6, 0.4, 0.5, 0.9);
    names(Y.trim) <- c("A","B","C","D","E","F","G","H");
    W <- rep(1,length(Y));
    W.trim <- rep(1,length(Y.trim));


    S.gpav  <- gpav(Y, W=NULL, adj); ## if W=NULL it is assumed that W is unitary
    S.gpav2 <- gpav(Y, W=W, adj);
    S.check <- c(0.53333333,0.30000000,0.53333333,0.53333333,0.53333333,0.53333333,0.53333333,0.55000000,0.55000000,0.90000000);
    names(S.check) <- c("H","I","J","F","D","B","C","E","G","A");

    expect_equal(S.gpav,  S.check, tolerance=1e-8);
    expect_equal(S.gpav2, S.check, tolerance=1e-8);
    expect_error(gpav(Y.trim,W=NULL,adj), "mismatch between the number of classes between Y and adj", fixed=TRUE);
    expect_error(gpav(Y,W=W.trim,adj), "mismatch between the number of classes between Y and adj", fixed=TRUE);
})

test_that("gpav.over.examples works", {
    S <- make.scores();
    g <- make.graph();

    S.gpav <- gpav.over.examples(S, g, W=NULL);

    pr1 <- c(0.90000000, 0.53333333, 0.53333333, 0.53333333, 0.55000000, 0.53333333, 0.55000000, 0.53333333, 0.30000000, 0.53333333);
    pr2 <- c(0.83333333, 0.83333333, 0.80000000, 0.83333333, 0.43333333, 0.63333333, 0.43333333, 0.43333333, 0.63333333, 0.63333333);
    pr3 <- c(0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.70000000, 0.10000000, 0.70000000);
    S.check <- rbind(pr1,pr2,pr3);
    colnames(S.check) <- LETTERS[1:length(pr1)];
    
    S.error <- S[,-which(colnames(S) %in% c("D","H"))];

    expect_equal(S.gpav, S.check, tolerance=1e-8);
    expect_error(gpav.over.examples(S.error, g, W=NULL), "mismatch between the number of nodes of the graph g and the number of classes of the scores matrix S");
})

test_that("gpav.parallel works", {
    S <- make.scores();
    g <- make.graph();

    if(Sys.info()['sysname']!="Windows"){ 
        S.gpav1core <- gpav.parallel(S, W=NULL, g, ncores=1);
        S.gpav2core <- gpav.parallel(S, W=NULL, g, ncores=2);
    }

    pr1 <- c(0.90000000, 0.53333333, 0.53333333, 0.53333333, 0.55000000, 0.53333333, 0.55000000, 0.53333333, 0.30000000, 0.53333333);
    pr2 <- c(0.83333333, 0.83333333, 0.80000000, 0.83333333, 0.43333333, 0.63333333, 0.43333333, 0.43333333, 0.63333333, 0.63333333);
    pr3 <- c(0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.70000000, 0.10000000, 0.70000000);
    S.check <- rbind(pr1,pr2,pr3);
    colnames(S.check) <- LETTERS[1:length(pr1)];
    
    S.error <- S[,-which(colnames(S) %in% c("D","H"))];

    expect_equal(S.gpav2core, S.check, tolerance=1e-8);
    expect_equal(S.gpav1core, S.check, tolerance=1e-8);
    expect_error(gpav.parallel(S.error, g, W=NULL), "mismatch between the number of nodes of the graph g and the number of classes of the scores matrix S");
})    

test_that("gpav.parallel works with ncores>2", {
    S <- make.scores();
    g <- make.graph();

    pr1 <- c(0.90000000, 0.53333333, 0.53333333, 0.53333333, 0.55000000, 0.53333333, 0.55000000, 0.53333333, 0.30000000, 0.53333333);
    pr2 <- c(0.83333333, 0.83333333, 0.80000000, 0.83333333, 0.43333333, 0.63333333, 0.43333333, 0.43333333, 0.63333333, 0.63333333);
    pr3 <- c(0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.70000000, 0.10000000, 0.70000000);
    S.check <- rbind(pr1,pr2,pr3);
    colnames(S.check) <- LETTERS[1:length(pr1)];

    ## apparently R CMD check only allows a maximum of two cores -> see a 'NB:...' in the R Packages Book (http://r-pkgs.had.co.nz/check.html)
    ## skip on CRAN -- avoid fail message 'simultaneous processes spawned'
    ## not_on_cran allows to execute test locally (with covr) but not on CRAN check
    not_on_cran <- function(){identical(Sys.getenv("NOT_CRAN"), "TRUE");}  ## set NOT_CRAN environment variable to TRUE
    if(not_on_cran()){
        S.gpav0core <- gpav.parallel(S, W=NULL, g, ncores=0);
        expect_equal(S.gpav0core, S.check, tolerance=1e-8); 
    }else{
        skip_on_cran(); ## skip test on CRAN
    }
})

# test_that("gpav.vanilla works", {
     ## to do
# })


# test_that("gpav.holdout works", {
    ## to do
# })








