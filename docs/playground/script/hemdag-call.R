#!/usr/bin/Rscript

## load library
library(HEMDAG);
suppressPackageStartupMessages(library(graph)); ## silence biocgenerics masked messages...
library(optparse);

## command line arguments
## for a detailed description, please see the manual: https://cran.r-project.org/web/packages/HEMDAG/HEMDAG.pdf
optionList <- list(
    make_option(c("-o", "--organism"), type="character", default="7227_drome",
        help="organism name in the form <taxon>_<name> (def. 7227_drome)"),
    make_option(c("-d", "--domain"), type="character", default="mf",
        help="go domain. It can be: bp, mf or cc (def. mf)"),
    make_option(c("-e", "--exptype"), type="character", default="ho",
        help="type of dataset on which run HEMDAG. It can be: ho (hold-out) or cv (cross-validated) -- def. ho"),
    make_option(c("-f", "--flat"), type="character", default="svmlinear",
        help="flat classifier"),
    make_option(c("-p", "--positive"), type="character", default="descendants",
        help="positive nodes selection. It can be: children or descendants. Skip this parameter if only topdown strategy is applied (def. descendants)"),
    make_option(c("-b", "--bottomup"), type="character", default="tau",
        help="bottomup strategy. It can be: none, threshold.free, threshold, weighted.threshold.free, weighted.threshold or tau. If none only topdown strategy is applied (def. tau)"),
    make_option(c("-t", "--topdown"), type="character", default="gpav",
        help="topdown strategy. It can be: htd or gpav (def. gpav)"),
    make_option(c("-c", "--threshold"), type="character", default="seq(from=0.1, to=0.9, by=0.1)",
        help="threshold for the choice of positive nodes. It can be a fixed value or an array of values (def. seq(from=0.1, to=0.9, by=0.1))"),
    make_option(c("-w", "--weight"), type="character", default="NULL",
        help="weight for the choice of positive nodes. It can be a fixed value or an array of values (def. NULL)"),
    make_option(c("-m", "--metric"), type="character", default="auprc",
        help="performance metric on which maximize the parametric ensemble algorithms. It can be: auprc or fmax (def. auprc)"),
    make_option(c("-r", "--round"), type="integer", default="3",
        help="number of rounding digits to be applied for choosing the best Fmax. To be used only if metric is set to fmax (def. 3)"),
    make_option(c("-s", "--seed"), type="integer", default="23",
        help="seed for the random generator to create folds (def. 23)"),
    make_option(c("-k", "--fold"), type="integer", default="5",
        help="number of folds for the cross validation (def. 5)"),
    make_option(c("-l", "--parallel"), type="logical", default=FALSE, action="store_true",
        help="should the sequential or parallel version of gpav be run? If flag -p is \"on\" the gpav parallel version is run. NB: only gpav can be run in parallel (def. FALSE)"),
    make_option(c("-n", "--cores"), type="integer", default="1",
        help="number of cores to use for the parallel execution of gpav (def. 1)"),
    make_option(c("-z", "--norm"), type="logical", default=FALSE, action="store_true",
        help="should the flat score matrix be normalized? If flag -p is \"on\" the input flat scores is normalized (def. FALSE)"),
    make_option(c("-y", "--normtype"), type="character", default="none",
        help="type of normalization. It can be maxnorm or qnorm (def. none)")
);

optParser <- OptionParser(option_list=optionList);
opt <- parse_args(optParser);

prefix    <- opt$organism;
organism  <- strsplit(prefix,"_")[[1]][2];
exptype   <- opt$exptype;
flat      <- opt$flat;
positive  <- opt$positive;
bottomup  <- opt$bottomup;
topdown   <- opt$topdown;
domain    <- opt$domain;
threshold <- eval(parse(text=opt$threshold));
weight    <- eval(parse(text=opt$weight));
metric    <- opt$metric;
round     <- opt$round;
seed      <- opt$seed;
kk        <- opt$fold;
parallel  <- opt$parallel;
cores     <- opt$cores;
norm      <- opt$norm;
normtype  <- opt$normtype;
if(normtype == "none")
    normtype <- NULL;

## hemdag algorithm to be displayed in output file name -> 18 iso/tpr-dag ensemble combinations + gpav + htd (tot 20 hemdag family)
if(positive=="children" && bottomup=="threshold.free" && topdown=="htd")
    hemdag.name <- "tprTF";
if(positive=="children" && bottomup=="threshold" && topdown=="htd")
    hemdag.name <- "tprT";
if(positive=="children" && bottomup=="weighted.threshold.free" && topdown=="htd")
    hemdag.name <- "tprW";
if(positive=="children" && bottomup=="weighted.threshold" && topdown=="htd")
    hemdag.name <- "tprwt";
if(positive=="descendants" && bottomup=="threshold.free" && topdown=="htd")
    hemdag.name <- "descensTF";
if(positive=="descendants" && bottomup=="threshold" && topdown=="htd")
    hemdag.name <- "descensT";
if(positive=="descendants" && bottomup=="weighted.threshold.free" && topdown=="htd")
    hemdag.name <- "descensW";
if(positive=="descendants" && bottomup=="weighted.threshold" && topdown=="htd")
    hemdag.name <- "descensWT";
if(positive=="descendants" && bottomup=="tau" && topdown=="htd")
    hemdag.name <- "descensTAU";
if(positive=="children" && bottomup=="threshold.free" && topdown=="gpav")
    hemdag.name <- "isotprTF";
if(positive=="children" && bottomup=="threshold" && topdown=="gpav")
    hemdag.name <- "isotprT";
if(positive=="children" && bottomup=="weighted.threshold.free" && topdown=="gpav")
    hemdag.name <- "isotprW";
if(positive=="children" && bottomup=="weighted.threshold" && topdown=="gpav")
    hemdag.name <- "isotprWT";
if(positive=="descendants" && bottomup=="threshold.free" && topdown=="gpav")
    hemdag.name <- "isodescensTF";
if(positive=="descendants" && bottomup=="threshold" && topdown=="gpav")
    hemdag.name <- "isodescensT";
if(positive=="descendants" && bottomup=="weighted.threshold.free" && topdown=="gpav")
    hemdag.name <- "isodescensW";
if(positive=="descendants" && bottomup=="weighted.threshold" && topdown=="gpav")
    hemdag.name <- "isodescensWT";
if(positive=="descendants" && bottomup=="tau" && topdown=="gpav")
    hemdag.name <- "isodescensTAU";
if(bottomup=="none" && topdown=="gpav")
    hemdag.name <- "gpav";
if(bottomup=="none" && topdown=="htd")
    hemdag.name <- "htd";

## I/O directories
data.dir <- paste0("../data/", exptype, "/");
res.dir  <- paste0("../res/",  exptype, "/");
if(!dir.exists(res.dir)){dir.create(res.dir, recursive=TRUE);}

## flat/ann/dag/testIndex files
files <- list.files(data.dir);
flat.file <- files[grep(paste0(organism, ".*", domain, ".scores.*", flat), files)];
ann.file  <- files[grep(paste0(domain,".ann"), files)];
dag.file  <- files[grep(paste0(domain,".dag"), files)];
if(exptype == "ho")
    idx.file  <- files[grep(paste0(domain,".testindex"), files)];

## load data
S <- get(load(paste0(data.dir, flat.file)));
g <- get(load(paste0(data.dir, dag.file)));
ann <- get(load(paste0(data.dir, ann.file)));
if(exptype == "ho")
    testIndex <- get(load(paste0(data.dir, idx.file)));

## shrink graph g to terms of matrix S -- if number of nodes between g and S mismatch
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

## elapsed time
start.elapsed <- proc.time();

## HEMDAG calling
if(exptype == "ho"){
    if(bottomup=="none"){
        if(topdown=="gpav"){
            S.hier <- gpav.holdout(S=S, g=g, testIndex=testIndex, W=NULL, parallel=parallel,
                ncores=cores, norm=norm, norm.type=normtype);
        }else{
            S.hier <- htd.holdout(S=S, g=g, testIndex=testIndex, norm=norm, norm.type=normtype);
        }
    }else{
        S.hier <- tpr.dag.holdout(S, g, ann=ann, testIndex=testIndex, norm=norm, norm.type=normtype,
            positive=positive, bottomup=bottomup, topdown=topdown, W=NULL, parallel=parallel, ncores=cores,
            threshold=threshold, weight=weight, kk=kk, seed=seed, metric=metric, n.round=round);
    }
}else{
    if(bottomup=="none"){
        if(topdown=="gpav"){
            S.hier <- gpav.vanilla(S=S, g=g, W=NULL, parallel=parallel, ncores=cores, norm=norm, norm.type=normtype);
        }else{
            S.hier <- htd.vanilla(S=S, g=g, norm=norm, norm.type=normtype);
        }
    }else{
        S.hier <- tpr.dag.cv(S, g, ann=ann, norm=norm, norm.type=normtype, positive=positive, bottomup=bottomup, 
            topdown=topdown, W=NULL, parallel=parallel, ncores=cores, threshold=threshold, weight=weight, 
            kk=kk, seed=seed, metric=metric, n.round=round);
    }
}

stop.elapsed <- proc.time() - start.elapsed;
timing.s <- stop.elapsed["elapsed"];
timing.m <- round(timing.s/(60),4);
timing.h <- round(timing.m/(60),4);
cat(hemdag.name, "running time:", timing.s["elapsed"], "(seconds)", "|", timing.m["elapsed"], "(minutes)", "|" , timing.h["elapsed"], "(hours)", "\n\n");

## store results
## outname
fname <- unlist(strsplit(flat.file, split="[.,_]"));
outname <- paste0(fname[-((length(fname)-1):length(fname))], collapse="_");
if(norm==TRUE && !(is.null(normtype)))
    outname <- paste0(outname,"_",normtype);
if(exptype == "ho"){
    save(S.hier, file=paste0(res.dir, outname, "_", hemdag.name, "_holdout.rda"), compress=TRUE);
}else{
    save(S.hier, file=paste0(res.dir, outname, "_", hemdag.name, ".rda"), compress=TRUE);
}

