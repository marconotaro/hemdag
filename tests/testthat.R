library(testthat);
library(HEMDAG);

test_check("HEMDAG");

## unneeded test. welcome message already tested in test_check
# test_that("start up message works", {
#     info <- sessionInfo();
#     if(is.null(info$otherPkgs$HEMDAG)){
#         expect_message(library(HEMDAG), "HEMDAG: Hierarchical Ensemble Methods for DAG-structured taxonomies\nPlease cite HEMDAG if you use it: see citation('HEMDAG') for details\n", fixed=TRUE);
#     }else{
#         detach("package:HEMDAG", unload=TRUE);
#         expect_message(library(HEMDAG), "HEMDAG: Hierarchical Ensemble Methods for DAG-structured taxonomies\nPlease cite HEMDAG if you use it: see citation('HEMDAG') for details\n", fixed=TRUE);
#     }
# })