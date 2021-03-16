#!/usr/bin/Rscript

start.time <- proc.time();

library(HEMDAG);
library(optparse);

optionList <- list(
    make_option(c("-o", "--organism"), type="character", default="7227_drome",
        help="organism name in the form <taxon>_<name> (def 7227_drome)"),
    make_option(c("-d", "--domain"),  type="character", default="mf",
        help="gene ontology domain -- bp, mf, cc (def. mf)"),
     make_option(c("-e", "--exptype"), type="character", default="cv",
        help="type of dataset on which run HEMDAG. It can be: ho (hold-out) or cv (cross-validated) -- def. ho"),
    make_option(c("-f", "--flat"), type="character", default="ranger",
        help="flat classifier (def. svmlinear)"),
    make_option(c("-a", "--algorithm"), type="character", default="isodescensTAU",
        help="hierarchical correction algorithm (def. isodescensTAU)")
);

optParser <- OptionParser(option_list=optionList);
opt <- parse_args(optParser);

## setting input argument
prefix    <- opt$organism;
taxon     <- strsplit(prefix,"_")[[1]][1];
organism  <- strsplit(prefix,"_")[[1]][2];
domain    <- opt$domain;
exptype   <- opt$exptype;
flat      <- opt$flat;
algorithm <- opt$algorithm;

## I/O
data.dir <- paste0("../data/", exptype, "/");
res.dir  <- paste0("../res/",  exptype, "/");
if(!dir.exists(res.dir)){dir.create(res.dir, recursive=TRUE);}

## flat/ann/dag/testIndex file
files <- list.files(data.dir);
flat.file <- files[grep(paste0(organism, ".*", domain, ".scores.*", flat), files)];
ann.file  <- files[grep(paste0(domain,".ann"), files)];
dag.file  <- files[grep(paste0(domain,".dag"), files)];
if(exptype == "ho")
    idx.file  <- files[grep(paste0(domain,".testindex"), files)];

## hier file
files <- list.files(res.dir);
hier.file <- files[grep(paste0(organism, ".*", domain, ".scores.*", flat, ".", algorithm), files)];

## load data
S <- get(load(paste0(data.dir, flat.file)));
S.hier <- get(load(paste0(res.dir, hier.file)));
g <- get(load(paste0(data.dir, dag.file)));
ann <- get(load(paste0(data.dir, ann.file)));
if(exptype == "ho")
    testIndex <- get(load(paste0(data.dir, idx.file)));

## if number of nodes between g and S mismatch -> shrink graph g to terms of matrix S
## eg. during flat learning you removed from S all those terms having less than N annotations 
root <- root.node(g);
nd <- colnames(S);
class.check <- ncol(S) != graph::numNodes(g);
if(class.check){
    root.check <- root %in% colnames(S);
    if(!root.check)
        nd <- c(root, nd);
    g <- build.subgraph(nd, g, edgemode="directed");
    ann <- ann[, colnames(S)];
}

## remove root node S nad S.hier score matrix (if any)
root <- root.node(g);
if((root %in% colnames(S)) && (root %in% colnames(S.hier))){
	S <- S[,-which(colnames(S)==root)];
    S.hier <- S.hier[,-which(colnames(S.hier)==root)];
	cat("root node removed from flat and hierarchical score matrix\n");
}

## remove root node from annotation matrix (if any)
if(root %in% colnames(ann)){
    ann <- ann[,-which(colnames(ann)==root)];
    cat("root node removed from annotation matrix\n");
}

## shrink S to testIndex
if(exptype == "ho"){
    S <- S[testIndex, ];
    ann <- ann[testIndex, colnames(S)];
}

## compute flat perf
auc.flat  <- auroc.single.over.classes(target=ann, predicted=S, folds=NULL, seed=NULL);
prc.flat  <- auprc.single.over.classes(target=ann, predicted=S, folds=NULL, seed=NULL);
fmax.flat <- compute.fmax(target=ann, predicted=S, n.round=3, b.per.example=TRUE, folds=NULL, seed=NULL, verbose=FALSE);
pxr.flat  <- precision.at.given.recall.levels.over.classes(target=ann, predicted=S, folds=NULL, seed=NULL, recall.levels=seq(from=0.1, to=1, by=0.1));
cat(taxon, organism, domain, flat, algorithm, "flat performance done\n");

## compute hier perf
auc.hier  <- auroc.single.over.classes(target=ann, predicted=S.hier, folds=NULL, seed=NULL);
prc.hier  <- auprc.single.over.classes(target=ann, predicted=S.hier, folds=NULL, seed=NULL);
fmax.hier <- compute.fmax(target=ann, predicted=S.hier, n.round=3, b.per.example=TRUE, folds=NULL, seed=NULL, verbose=FALSE);
pxr.hier  <- precision.at.given.recall.levels.over.classes(target=ann, predicted=S.hier, folds=NULL, seed=NULL, recall.levels=seq(from=0.1, to=1, by=0.1));
cat(taxon, organism, domain, flat, algorithm, "hierarchical performance done\n");

## storing
outname <- gsub("scores", "perfmeas", hier.file)
save(auc.flat, prc.flat, fmax.flat, pxr.flat, auc.hier, prc.hier, fmax.hier, pxr.hier, file=paste0(res.dir, outname), compress=TRUE);

## timing
timing.s <- proc.time() - start.time;
timing.m <- round(timing.s/(60),4);
timing.h <- round(timing.m/(60),4);
cat("\n");
cat("elapsed time:", timing.s["elapsed"], "(seconds)", "|", timing.m["elapsed"], "(minutes)", "|" , timing.h["elapsed"], "(hours)", "\n\n");

