library(HEMDAG);

context("test normalization function");

make.scores <- function(){
    pr1 <- c(0.9, 0.4, 0.4, 0.3, 0.6, 0.4, 0.5, 0.9, 0.3, 0.8);
    pr2 <- c(0.7, 0.9, 0.8, 0.9, 0.5, 0.2, 0.3, 0.5, 0.7, 1.0);
    pr3 <- c(0.1, 0.9, 0.7, 1.0, 0.9, 0.4, 0.8, 0.9, 0.1, 0.8);
    S <- rbind(pr1,pr2,pr3);
    colnames(S) <- LETTERS[1:length(pr1)];
    return(S);
}

test_that("Do.flat.scores.normalization works",{
    S <- make.scores();
    tmpdir <- paste0(tempdir(),"/");
    save(S, file=paste0(tmpdir,"scores.rda"));
    flat.dir <- flat.norm.dir <- tmpdir;

    Do.flat.scores.normalization(norm.type="MaxNorm", flat.file="scores", flat.dir=flat.dir, flat.norm.dir=flat.norm.dir);
    Do.flat.scores.normalization(norm.type="Qnorm", flat.file="scores", flat.dir=flat.dir, flat.norm.dir=flat.norm.dir);

    S.maxnorm.check <- scores.normalization(norm.type="MaxNorm", S);
    S.qnorm.check <- scores.normalization(norm.type="Qnorm", S);

    S.maxnorm <- get(load(paste0(tmpdir,"MaxNorm.scores.rda")));
    S.qnorm <- get(load(paste0(tmpdir,"Qnorm.scores.rda")));

    expect_equal(S.maxnorm, S.maxnorm.check);
    expect_equal(S.qnorm, S.qnorm.check);
})
