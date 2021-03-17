#!/usr/bin/Rscript

start.time <- proc.time();

library(HEMDAG);
library(optparse);

optionList <- list(
    make_option(c("-e", "--exptype"), type="character", default="ho",
        help="type of dataset on which run HEMDAG. It can be: ho (hold-out) or cv (cross-validated) -- def. ho")
);

optParser <- OptionParser(option_list=optionList);
opt <- parse_args(optParser);

exptype   <- opt$exptype;

data.dir <- paste0("../data/", exptype, "/");
res.dir  <- paste0("../res/",  exptype, "/");
data.files <- list.files(data.dir);
res.files  <- list.files(res.dir);

orgs  <- c("7227_drome"); ## organism(s) list
flats <- c("svmlinear");  ## flat classifier(s) list
algs  <- c("gpav", "isotprTF", "isotprW", "isodescensTF", "isodescensW", "isodescensTAU");  ## HEMDAG algorithm(s) list
onts  <- c("mf");   ## GO domain(s) list (bp, mf, cc)

for(org in orgs){
    for(flat in flats){
        for(alg in algs){
            for(ont in onts){
                dag.file  <- data.files[grep(paste0(ont,".dag"), data.files)];
                hier.file <- res.files[grep(paste0(org, ".*", ont, ".scores.*", flat, ".", alg), res.files)];
                g <- get(load(paste0(data.dir, dag.file)));
                S.hier <- get(load(paste0(res.dir, hier.file)));
                root <- root.node(g);
                g <- build.subgraph(colnames(S.hier), g);
                check <- check.hierarchy(S.hier, g, root);
                if(check$status=="OK"){
                    cat(org, ont, paste0(flat,'+',alg), "check passed :)","\n");
                }else{
                    cat(org, ont, paste0(flat,'+',alg), "check failed :(","\n");
                }
            }
        }
    }
}

timing.s <- proc.time() - start.time;
timing.m <- round(timing.s/(60),4);
timing.h <- round(timing.m/(60),4);
cat("\n");
cat("elapsed time:", timing.s["elapsed"], "(seconds)", "|", timing.m["elapsed"], "(minutes)", "|" , timing.h["elapsed"], "(hours)", "\n\n");

