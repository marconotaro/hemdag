#!/usr/bin/Rscript

## script to compute auprc / auroc / pxr / fmax one-shot --- flat vs hier

## start perf
start.time <- proc.time();

library(HEMDAG);
library(optparse);

optionList <- list(
    make_option(c("-o", "--organism"), type="character", default="7227.drome",
        help="organism name in the form <taxon>.<name>"),
    make_option(c("-d", "--domain"),  type="character", default="mf",
        help="gene ontology domain -- bp, mf, cc"),
    make_option(c("-f", "--flat"), type="character", default="svmlinear",
        help="flat classifier"),
    make_option(c("-a", "--algorithm"), type="character", default="isodescensTAU",
        help="hierarchical algorithm")
);

optParser <- OptionParser(option_list=optionList);
opt <- parse_args(optParser);

## setting input argument
org <- opt$organism;
ont <- opt$domain;
alg <- opt$flat;
met <- opt$algorithm;

## I/O
data.dir <- paste0("../data/ho/");
res.dir <- "../res/ho/";

g         <- get(load(paste0(data.dir, org,".go.",ont,".dag.20dec17-16jun20.rda")));
testIndex <- get(load(paste0(data.dir, org,".go.",ont,".testindex.20dec17-16jun20.rda")));
ann       <- get(load(paste0(data.dir, org,".go.",ont,".ann.20dec17-16jun20.rda")));
S         <- get(load(paste0(data.dir, org,".go.",ont,".scores.",alg,".holdout.rda")));
S.hier    <- get(load(paste0(res.dir,  org,".go.",ont,".scores.",alg,".",met,".holdout.rda")));

## remove root node from S and S.hier scores (if any)
root <- root.node(g);
if(root %in% colnames(S)){
	S <- S[,-which(colnames(S)==root)];
	cat("root node removed from flat scores matrix S\n");
}
if(root %in% colnames(S.hier)){
	S.hier <- S.hier[,-which(colnames(S.hier)==root)];
	cat("root node removed from hierarchical scores matrix S.hier\n");
}

## shrink S to testIndex
S <- S[testIndex, ];
ann <- ann[testIndex, colnames(S)];

## compute flat perf
auc.flat  <- auroc.single.over.classes(target=ann, predicted=S, folds=NULL, seed=NULL);
prc.flat  <- auprc.single.over.classes(target=ann, predicted=S, folds=NULL, seed=NULL);
fmax.flat <- compute.fmax(target=ann, predicted=S, n.round=3, b.per.example=TRUE, folds=NULL, seed=NULL, verbose=FALSE);
pxr.flat  <- precision.at.given.recall.levels.over.classes(target=ann, predicted=S, folds=NULL, seed=NULL, recall.levels=seq(from=0.1, to=1, by=0.1));
cat(tax, org, ont, alg, met, "flat performance done\n");

## compute hier perf
auc.hier  <- auroc.single.over.classes(target=ann, predicted=S.hier, folds=NULL, seed=NULL);
prc.hier  <- auprc.single.over.classes(target=ann, predicted=S.hier, folds=NULL, seed=NULL);
fmax.hier <- compute.fmax(target=ann, predicted=S.hier, n.round=3, b.per.example=TRUE, folds=NULL, seed=NULL, verbose=FALSE);
pxr.hier  <- precision.at.given.recall.levels.over.classes(target=ann, predicted=S.hier, folds=NULL, seed=NULL, recall.levels=seq(from=0.1, to=1, by=0.1));
cat(org, ont, alg, met, "hierarchical performance done\n");

## storing
save(auc.flat, prc.flat, fmax.flat, pxr.flat, auc.hier, prc.hier, fmax.hier, pxr.hier, file=paste0(res.dir, tax,".",org,".go.",ont,".perfmeas.",alg,".",met,".holdout.rda"), compress=TRUE);

## timing
timing.s <- proc.time() - start.time;
timing.m <- round(timing.s/(60),4);
timing.h <- round(timing.m/(60),4);
cat("\n");
cat("elapsed time:", timing.s["elapsed"], "(seconds)", "|", timing.m["elapsed"], "(minutes)", "|" , timing.h["elapsed"], "(hours)", "\n\n");

