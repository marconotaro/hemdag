library(HEMDAG);

context("test annotation utility functions");

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

make.adj <- function(){
    pr1 <- c(0.73,0.87,0.93,0.77);
    pr2 <- c(0.89,0.98,0.74,0.95);
    pr4 <- c(0.71,0.92,0.83,0.99);
    pr5 <- c(0.75,0.64,0.78,0.88);
    W <- rbind(pr1,pr2,pr4,pr5);
    colnames(W) <- c("pr1","pr2","pr4","pr5");
    return(W);
}

make.notca.ann <- function(){
    pr1 <- c(1,1,1,1,0,1,0,0,1,1);
    pr2 <- c(1,1,1,1,1,1,1,1,0,0);
    pr3 <- c(1,1,1,1,0,1,0,0,1,0);
    ann <- rbind(pr1,pr2,pr3);
    colnames(ann) <- c("A","B","C","D","E","F","G","H","I","J");
    return(ann);
} 

make.spec.ann <- function(){
    pr1 <- c(1,0,1);
    pr2 <- c(0,1,0);
    pr3 <- c(1,0,1);
    spec.ann <- rbind(pr1,pr2,pr3);
    colnames(spec.ann) <- c("C","H","J");
    return(spec.ann);
}

test_that("specific.annotation.matrix works", {
    ann.prname <- "pr1 A|C\npr2 B";
    ann.entrex <- "1 A|C\n2 B";

    file.txt <- paste0(tempdir(),"/","annfile.txt");
    file.zip <- paste0(tempdir(),"/","annfile.txt.gz");
    file.etx <- paste0(tempdir(),"/","annfile.entrez");

    writeLines(ann.prname, con=file.txt);
    writeLines(ann.prname, con=file.zip);
    writeLines(ann.entrex, con=file.etx);

    pr1 <- c(1,0,1);
    pr2 <- c(0,1,0);
    ann.goldID <- ann.gold <- rbind(pr1,pr2);
    colnames(ann.gold) <- c("A","B","C");
    rownames(ann.goldID) <- 1:nrow(ann.gold);
    colnames(ann.goldID) <- colnames(ann.gold);

    ann.txt <- specific.annotation.matrix(file=file.txt); 
    ann.zip <- specific.annotation.matrix(file=file.zip); 
    ann.etx <- specific.annotation.matrix(file=file.etx); 

    expect_equal(ann.txt, ann.gold);
    expect_equal(ann.zip, ann.gold);
    expect_equal(ann.etx, ann.goldID);
})

test_that("specific.annotation.list works", {
    spec.ann <- make.spec.ann();
    ann.list <- specific.annotation.list(spec.ann);

    expect_equal(ann.list[["pr1"]], c("C","J"));
    expect_equal(ann.list[["pr2"]], "H");
    expect_equal(ann.list[["pr3"]], c("C","J"));
})

test_that("transitive.closure.annotations works", {
    g <- make.graph();
    anc <- build.ancestors(g);
    spec.ann <- make.spec.ann();
    tca <- transitive.closure.annotations(spec.ann, anc);
    
    pr1 <- c(1,1,1,1,0,1,0,0,1);
    pr2 <- c(1,1,1,1,1,1,1,1,0);
    pr3  <- c(1,1,1,1,0,1,0,0,1);
    tca.control <- rbind(pr1,pr2,pr3);
    colnames(tca.control) <- c("A","B","C","D","E","F","G","H","J");

    expect_equal(tca, tca.control);
})

test_that("full.annotation.matrix works", {
    g <- make.graph();
    anc <- build.ancestors(g);
    W <- make.adj();
    spec.ann <- make.spec.ann();
    ann <- make.notca.ann();

    full.ann <- full.annotation.matrix(W,anc,spec.ann); ## tca performed

    pr1 <- c(1,1,1,1,0,1,0,0,1);
    pr2 <- c(1,1,1,1,1,1,1,1,0);
    pr4 <- c(0,0,0,0,0,0,0,0,0);
    pr5 <- c(0,0,0,0,0,0,0,0,0);
    full.ann.control <- rbind(pr1,pr2,pr4,pr5);
    colnames(full.ann.control) <- c("A","B","C","D","E","F","G","H","J");

    expect_equal(full.ann, full.ann.control);
})

test_that("build.submatrix works", {
    ann <- make.notca.ann();
    subann <- build.submatrix(ann,2);

    pr1 <- c(1,1,1,1,1);
    pr2 <- c(1,1,1,1,1);
    pr3 <- c(1,1,1,1,1);
    control <- rbind(pr1,pr2,pr3);
    colnames(control) <- c("A","B","C","D","F");

    expect_equal(subann, control);
})

test_that("check.annotation.matrix.integrity works", {
    g <- make.graph();
    anc <- build.ancestors(g);
    spec.ann <- make.spec.ann();
    ann <- make.notca.ann();
    tca <- transitive.closure.annotations(spec.ann, anc);

    expect_output(check.annotation.matrix.integrity(anc, spec.ann, ann), "check.annotation.matrix: NOT_OK. Transitive closure NOT RESPECTED");
    expect_output(check.annotation.matrix.integrity(anc, spec.ann, tca), "check.annotation.matrix: OK ");
})


