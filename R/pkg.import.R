#' @import graph
#' @import RBGL
#' @import  precrec				   
#' @import  PerfMeas			   
#' @import  preprocessCore
#' @import methods
#' @import iterators
#' @import parallel
#' @import foreach
#' @import doParallel
#' @importFrom plyr mapvalues
#' @importFrom utils read.table write.table
#' @useDynLib HEMDAG, .registration=TRUE 


## Quiet concerns of R CMD check. 
## Avoid this warning: no visible binding for global variable
if(getRversion() >= "2.15.1"){
	utils::globalVariables(c("curvetypes","detectCores", "i"));
}
