#!/usr/bin/Rscript

suppressPackageStartupMessages(library(HEMDAG));
library(optparse);

optionList <- list(
    make_option(c("-o", "--organism"), type="character", default="7227_drome",
        help="organism name in the form <taxon>_<name> (def 7227_drome)"),
    make_option(c("-d", "--domain"),  type="character", default="mf",
        help="gene ontology domain -- bp, mf, cc (def. mf)"),
     make_option(c("-e", "--exptype"), type="character", default="ho",
        help="type of dataset on which run HEMDAG. It can be: ho (hold-out) or cv (cross-validated) -- def. ho"),
    make_option(c("-f", "--flat"), type="character", default="svmlinear",
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
data.files <- list.files(data.dir);
res.files  <- list.files(res.dir);
if(!dir.exists(res.dir)){dir.create(res.dir, recursive=TRUE);}

## load data
dag.file  <- data.files[grep(paste0(domain,".dag"), data.files)];
hier.file <- res.files[grep(paste0(organism, ".*", domain, ".scores.*", flat, ".", algorithm), res.files)];

## check if flat|ann|dag exists
if(length(dag.file)==0 || length(hier.file)==0)
    stop("no dag|hier file found\n");

## check constraints violation
g <- get(load(paste0(data.dir, dag.file)));
S.hier <- get(load(paste0(res.dir, hier.file)));
root <- root.node(g);
g <- build.subgraph(colnames(S.hier), g);
check <- check.hierarchy(S.hier, g, root);
if(check$status=="OK"){
    cat(taxon, organism, domain, paste0(flat,'+',algorithm), "check passed :)","\n");
}else{
    cat(taxon, organism, domain, paste0(flat,'+',algorithm), "check failed :(","\n");
}
