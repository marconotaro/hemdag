library(HEMDAG);

context("test GPAV");

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

test_that("GPAV works", {
    g <- make.graph();
    adj <- adj.upper.tri(g);

    Y <- c(0.9, 0.4, 0.4, 0.3, 0.6, 0.4, 0.5, 0.9, 0.3, 0.8);
    names(Y) <- c("A","B","C","D","E","F","G","H","I","J");
    Y.trim <- c(0.9, 0.4, 0.4, 0.3, 0.6, 0.4, 0.5, 0.9);
    names(Y.trim) <- c("A","B","C","D","E","F","G","H");
    W <- rep(1,length(Y));
    W.trim <- rep(1,length(Y.trim));


    S.gpav  <- GPAV(Y, W=NULL, adj); ## if W=NULL it is assumed that W is unitary
    S.gpav2 <- GPAV(Y, W=W, adj);
    S.check <- c(0.53333333,0.30000000,0.53333333,0.53333333,0.53333333,0.53333333,0.53333333,0.55000000,0.55000000,0.90000000);
    names(S.check) <- c("H","I","J","F","D","B","C","E","G","A");

    expect_equal(S.gpav,  S.check, tolerance=1e-8);
    expect_equal(S.gpav2, S.check, tolerance=1e-8);
    expect_error(GPAV(Y.trim,W=NULL,adj), "GPAV: mismatch between the number of classes between 'Y' and 'adj'", fixed=TRUE);
    expect_error(GPAV(Y,W=W.trim,adj), "GPAV: mismatch between the number of classes between 'Y' and 'adj'", fixed=TRUE);
})

test_that("GPAV.over.examples works", {
    S <- make.scores();
    g <- make.graph();

    S.gpav <- GPAV.over.examples(S, g, W=NULL);

    pr1 <- c(0.90000000, 0.53333333, 0.53333333, 0.53333333, 0.55000000, 0.53333333, 0.55000000, 0.53333333, 0.30000000, 0.53333333);
    pr2 <- c(0.83333333, 0.83333333, 0.80000000, 0.83333333, 0.43333333, 0.63333333, 0.43333333, 0.43333333, 0.63333333, 0.63333333);
    pr3 <- c(0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.70000000, 0.10000000, 0.70000000);
    S.check <- rbind(pr1,pr2,pr3);
    colnames(S.check) <- LETTERS[1:length(pr1)];
    
    S.error <- S[,-which(colnames(S) %in% c("D","H"))];

    expect_equal(S.gpav, S.check, tolerance=1e-8);
    expect_error(GPAV.over.examples(S.error, g, W=NULL), 
        "GPAV: mismatch between the number of nodes of the graph and the number of classes of the scores matrix");
})

test_that("GPAV.parallel works", {
    S <- make.scores();
    g <- make.graph();

    if(Sys.info()['sysname']!="Windows"){ 
        S.gpav1core <- GPAV.parallel(S, W=NULL, g, ncores=1);
        S.gpav2core <- GPAV.parallel(S, W=NULL, g, ncores=2); 
    }

    pr1 <- c(0.90000000, 0.53333333, 0.53333333, 0.53333333, 0.55000000, 0.53333333, 0.55000000, 0.53333333, 0.30000000, 0.53333333);
    pr2 <- c(0.83333333, 0.83333333, 0.80000000, 0.83333333, 0.43333333, 0.63333333, 0.43333333, 0.43333333, 0.63333333, 0.63333333);
    pr3 <- c(0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.74000000, 0.70000000, 0.74000000, 0.70000000, 0.10000000, 0.70000000);
    S.check <- rbind(pr1,pr2,pr3);
    colnames(S.check) <- LETTERS[1:length(pr1)];
    
    S.error <- S[,-which(colnames(S) %in% c("D","H"))];

    expect_equal(S.gpav2core, S.check, tolerance=1e-8);
    expect_equal(S.gpav1core, S.check, tolerance=1e-8);
    expect_error(GPAV.parallel(S.error, g, W=NULL),
        "GPAV: mismatch between the number of nodes of the graph and the number of classes of the scores matrix");
})    

test_that("GPAV.parallel works with ncores>2", {
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
        S.gpav0core <- GPAV.parallel(S, W=NULL, g, ncores=0);
        expect_equal(S.gpav0core, S.check, tolerance=1e-8); 
    }else{
        skip_on_cran(); ## skip test on CRAN
    }
})

test_that("Do.GPAV works", {
    g <- make.graph();
    S <- make.scores();
    ann <- make.ann();
    S.noroot <- S[,-which(colnames(S)==root.node(g))];
    S.shrank <- S[,colnames(do.submatrix(ann,2))];
    
    tmpdir <- paste0(tempdir(),"/");
    save(g, file=paste0(tmpdir,"graph.rda"));
    save(ann, file=paste0(tmpdir,"labels.rda"));
    save(S, file=paste0(tmpdir,"scores.rda"));
    save(S.noroot, file=paste0(tmpdir,"scores-noroot.rda"));
    save(S.shrank, file=paste0(tmpdir,"scores-shrunk.rda"));
    dag.dir <- flat.dir <- ann.dir <- tmpdir;
    hierScore.dir <- perf.dir <- tmpdir;
    recall.levels <- seq(from=0.25, to=1, by=0.25);
    dag.file <- "graph";
    flat.file <- "scores";
    flat.file.noroot <- "scores-noroot";
    flat.file.shrunk <- "scores-shrunk";
    ann.file <- "labels";
    
    ## test S -> NoNorm
    expect_output(
        Do.GPAV(norm=TRUE, norm.type=NULL, W=NULL, parallel=FALSE, ncores=1, folds=NULL, seed=23, n.round=3, 
            f.criterion ="F", recall.levels=recall.levels, compute.performance=TRUE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir), 
        "FLAT PERFORMANCE: DONE\\nHIERARCHICAL CORRECTION: DONE\\nHIERARCHICAL PERFORMANCE: DONE"
    );

    ## test hierarchical scores matrix
    S.hier <- get(load(paste0(tmpdir,"scores.hierScores.GPAV.rda")));
    S.gpav <- GPAV.over.examples(S, g, W=NULL); 
    expect_equal(S.hier, S.gpav);

    ## testing flat vs hierarchical performance
    perfs <- mget(load(paste0(tmpdir,"PerfMeas.scores.hierScores.GPAV.rda")));
    ann.noroot <- ann[,-which(colnames(S)==root.node(g))];
    S.gpav.noroot <- S.gpav[,-which(colnames(S.gpav)==root.node(g))];

    AUC.flat <- AUROC.single.over.classes(ann.noroot, S.noroot, folds=NULL, seed=23);
    PRC.flat <- AUPRC.single.over.classes(ann.noroot, S.noroot, folds=NULL, seed=23);
    FMM.flat <- compute.Fmeasure.multilabel(ann.noroot, S.noroot, n.round=3, f.criterion="F", verbose=FALSE, b.per.example=TRUE, folds=NULL, seed=23);
    AUC.gpav <- AUROC.single.over.classes(ann.noroot, S.gpav.noroot, folds=NULL, seed=23);
    PRC.gpav <- AUPRC.single.over.classes(ann.noroot, S.gpav.noroot, folds=NULL, seed=23);
    FMM.gpav <- compute.Fmeasure.multilabel(ann.noroot, S.gpav.noroot, n.round=3, f.criterion="F", verbose=FALSE, b.per.example=TRUE, folds=NULL, seed=23);

    expect_equal(perfs$AUC.flat, AUC.flat);
    expect_equal(perfs$PRC.flat, PRC.flat);
    expect_equal(perfs$FMM.flat, FMM.flat);
    expect_equal(perfs$AUC.hier, AUC.gpav);
    expect_equal(perfs$PRC.hier, PRC.gpav);
    expect_equal(perfs$FMM.hier, FMM.gpav);

    ## test S -> MaxNorm 
    expect_output( 
        Do.GPAV(norm=FALSE, norm.type="MaxNorm", W=NULL, parallel=FALSE, ncores=1, folds=NULL, seed=23, n.round=3, 
            f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL), 
        "MaxNorm NORMALIZATION: DONE\\nHIERARCHICAL CORRECTION: DONE"
    );

    ## test S -> Qnorm
    expect_output(
        Do.GPAV(norm=FALSE, norm.type="Qnorm", W=NULL, parallel=FALSE, ncores=1, folds=NULL, seed=23, n.round=3, 
            f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL), 
        "Qnorm NORMALIZATION: DONE\\nHIERARCHICAL CORRECTION: DONE"
    );

    ## test misspelled norm method
    expect_error(
        Do.GPAV(norm=FALSE, norm.type="BlaBla", W=NULL, parallel=TRUE, ncores=2, folds=NULL, seed=23, n.round=3,
            f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL),
        "scores.normalization: the chosen normalization method is not among those available or it has been misspelled"
    );

    ## test misspelled performance
    expect_error(
        Do.GPAV(norm=FALSE, norm.type="Qnorm", W=NULL, parallel=TRUE, ncores=2, folds=NULL, seed=23, n.round=3,
            f.criterion="Fmax", recall.levels=recall.levels, compute.performance=TRUE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir),
        "GPAV: value of parameter 'f.criterion' misspelled"
    );

    ## test error norm=FALSE and norm.type=NULL
    expect_error(
        Do.GPAV(norm=FALSE, norm.type=NULL, W=NULL, parallel=TRUE, ncores=2, folds=NULL, seed=23, n.round=3,
            f.criterion="F", recall.levels=recall.levels, compute.performance=TRUE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir),
        "GPAV: If norm is set to FALSE, you need to specify a normalization method among those available"
    );

    ## test S without root 
    expect_output(
        Do.GPAV(norm=TRUE, norm.type=NULL, W=NULL, parallel=FALSE, ncores=1, folds=NULL, seed=23, n.round=3, 
            f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file.noroot, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL), 
        "HIERARCHICAL CORRECTION: DONE"
    );

    ## test GPAV parallel
    expect_output(
        Do.GPAV(norm=TRUE, norm.type=NULL, W=NULL, parallel=TRUE, ncores=2, folds=NULL, seed=23, n.round=3, 
            f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL), 
        "HIERARCHICAL CORRECTION: DONE"
    );

    ## test GPAV parallel warning 
    expect_output( 
        expect_warning(
            Do.GPAV(norm=TRUE, norm.type=NULL, W=NULL, parallel=TRUE, ncores=1, folds=NULL, seed=23, n.round=3,
                f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
                dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL),  
                "GPAV: set ncores greater than 2 to exploit the GPAV parallel version"
        ),
        "HIERARCHICAL CORRECTION: DONE"
    );

    expect_output( 
        expect_warning(
            Do.GPAV(norm=TRUE, norm.type=NULL, W=NULL, parallel=FALSE, ncores=2, folds=NULL, seed=23, n.round=3,
                f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
                dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL),  
                "GPAV: no GPAV parallel version is running, but ncores is higher or equal to 2.", " Set 'ncores' to 1 to run the sequential version or set 'parallel' to TRUE to run the parallel version"
        ),
        "HIERARCHICAL CORRECTION: DONE"
    );

    ## test normalization warning
    expect_output( 
        expect_warning(
            Do.GPAV(norm=TRUE, norm.type="MaxNorm", W=NULL, parallel=TRUE, ncores=2, folds=NULL, seed=23, n.round=3,
                f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file, ann.file=ann.file, 
                dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL),  
                "GPAV: If norm is set to TRUE, the input flat matrix is already normalized. Set norm.type to NULL and not to 'MaxNorm' to avoid this warning message"
        ),
        "HIERARCHICAL CORRECTION: DONE"
    );

    ## test compute.performance=TRUE and norm=FALSE
    expect_output(
        Do.GPAV(norm=FALSE, norm.type="MaxNorm", W=NULL, parallel=TRUE, ncores=2, folds=NULL, seed=23, n.round=3, 
            f.criterion="F", recall.levels=recall.levels, compute.performance=TRUE, flat.file=flat.file, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir), 
        "HIERARCHICAL CORRECTION: DONE"
    );

    ## test shrunk flat scores matrix -> matrix having a subset of ontology terms (on the basis of annotations number)
    expect_output(
        Do.GPAV(norm=TRUE, norm.type=NULL, W=NULL, parallel=TRUE, ncores=2, folds=NULL, seed=23, n.round=3, 
            f.criterion=NULL, recall.levels=NULL, compute.performance=FALSE, flat.file=flat.file.shrunk, ann.file=ann.file, 
            dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=NULL), 
        "HIERARCHICAL CORRECTION: DONE"
    );
})


# test_that("Do.GPAV.holdout works", {
    ## to be added
# })








