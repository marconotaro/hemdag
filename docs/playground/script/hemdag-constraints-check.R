#!/usr/bin/Rscript

## script to check scores consistency
library(HEMDAG);

data.dir <- "../data/ho/";
res.dir <- "../res/ho/";

orgs <- c("7227.drome");
algs <- c("svmlinear");
mets <- c("gpav", "isotprTF", "isotprW", "isodescensTF", "isodescensW", "isodescensTAU");
onts <- c("mf");

for(org in orgs){
    tax <- strsplit(org, split="[.]")[[1]][1];
    org <- strsplit(org, split="[.]")[[1]][2];
    for(alg in algs){
        for(met in mets){
            for(ont in onts){
                g <- get(load(paste0(data.dir, tax,".",org,".go.",ont,".dag.20dec17-16jun20.rda")));
                S.hier <- get(load(paste0(res.dir, tax,".",org,".go.",ont,".scores.",alg,".",met,".holdout.rda")));
                root <- root.node(g);
                g <- build.subgraph(colnames(S.hier), g);
                check <- check.hierarchy(S.hier, g, root);
                if(check$status=="OK"){
                    cat(org, ont, paste0(alg,'+',met), "check passed :)","\n");
                }else{
                    cat(org, ont, paste0(alg,'+',met), "check failed :(","\n");
                }
            }
        }
    }
}

