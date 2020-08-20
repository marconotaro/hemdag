library(HEMDAG);

context("test graph utility functions");

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
    ## optional: plotting the dag g
    ## library(Rgraphviz); plot(g);
}

test_that("root.node works", {
    g <- make.graph();
    root <- root.node(g);

    expect_equal(root, "A");
    expect_equal(length(root), 1);
})

test_that("compute.flipped.graph works", {
    if (!check.graph()){ skip("graph package cannot be loaded"); }

    g <- make.graph();
    ndg <- graph::numNodes(g);
    edg <- graph::numEdges(g);

    og <- compute.flipped.graph(g);
    ndog <- numNodes(og);
    edog <- numEdges(og);

    expect_s4_class(og,"graphNEL");
    expect_equal(ndg, 10);
    expect_equal(edg, 16);
    expect_equal(ndog, 10);
    expect_equal(edog, 16);
})

test_that("graph.levels works", {
    g <- make.graph();
    root <- root.node(g);
    lev <- graph.levels(g, root=root);

    expect_equal(lev[["level_0"]], "A");
    expect_equal(lev[["level_1"]], c("B","C","G"));
    expect_equal(lev[["level_2"]], c("D","E"));
    expect_equal(lev[["level_3"]], "F");
    expect_equal(lev[["level_4"]], c("H","I","J"));
    expect_equal(length(lev), 5);

    expect_error(graph.levels(g, root="R"), "root node not found in g. Insert the root node", fixed=TRUE);
    expect_error(graph.levels(g, root="B"), "root is not the right root node of g. Use the function root.node(g) to find the root node of g", fixed=TRUE);
})

test_that("get.parents works", {
    g <- make.graph();
    nd <- nodes(g);
    root <- root.node(g)
    parents <- get.parents(g, root=root);
    lev <- graph.levels(g, root=root);
    parents.tod <- get.parents.top.down(g, lev, root=root);
    parents.bup <- get.parents.bottom.up(g, lev, root=root);
    parents.tsort <- get.parents.topological.sorting(g, root=root);

    expect_equal(length(parents), 9);
    expect_equal(length(parents.tod), 9);
    expect_equal(length(parents.bup), 9);
    expect_equal(names(parents), c("B","C","D","E","F","G","H","I","J")); 
    expect_equal(names(parents.tod), c("B","C","G","D","E","F","H","I","J")); 
    expect_equal(names(parents.bup), c("J","I","H","F","E","D","G","C","B"));
    expect_equal(names(parents.tsort), c("G","E","C","B","D","F","J","I","H"));
    expect_equal(parents[["B"]], "A");
    expect_equal(parents[["C"]], "A");
    expect_equal(parents[["D"]], "B");
    expect_equal(parents[["E"]], c("A","G"));
    expect_equal(parents[["F"]], c("B","C","D"));
    expect_equal(parents[["G"]], "A");
    expect_equal(parents[["H"]], c("C","E","F","G"));
    expect_equal(parents[["I"]], c("D","F"));
    expect_equal(parents[["J"]], "F");

    expect_error(get.parents(g, root="R"), "root node not found in g. Insert the root node", fixed=TRUE);
    expect_error(get.parents(g, root="B"), "root is not the right root node of g. Use the function root.node(g) to find the root node of g", fixed=TRUE);
    expect_error(get.parents.top.down(g, root="R"), "root node not found in g. Insert the root node", fixed=TRUE);
    expect_error(get.parents.top.down(g, root="B"), "root is not the right root node of g. Use the function root.node(g) to find the root node of g", fixed=TRUE);
    expect_error(get.parents.bottom.up(g,root="R"), "root node not found in g. Insert the root node", fixed=TRUE);
    expect_error(get.parents.bottom.up(g,root="B"), "root is not the right root node of g. Use the function root.node(g) to find the root node of g", fixed=TRUE);
    expect_error(get.parents.topological.sorting(g,root="R"), "root node not found in g. Insert the root node", fixed=TRUE);
    expect_error(get.parents.topological.sorting(g,root="B"), "root is not the right root node of g. Use the function root.node(g) to find the root node of g",
        fixed=TRUE);
})

test_that("get.children works", {
    g <- make.graph();
    nd <- nodes(g);
    root <- root.node(g)
    children <- build.children(g);
    lev <- graph.levels(g, root=root);
    children.tod <- get.children.top.down(g,lev);
    children.bup <- get.children.bottom.up(g,lev);

    expect_equal(length(children), 10);
    expect_equal(length(children.tod), 10);
    expect_equal(length(children.bup), 10);
    expect_equal(names(children), c("A","B","C","D","E","F","G","H","I","J")); 
    expect_equal(names(children.tod), c("A","B","C","G","D","E","F","H","I","J")); 
    expect_equal(names(children.bup), c("J","I","H","F","E","D","G","C","B","A"));

    expect_equal(children[["A"]], c("B","C","G","E"));
    expect_equal(children[["B"]], c("D","F"));
    expect_equal(children[["C"]], c("F","H"));
    expect_equal(children[["D"]], c("F","I"));
    expect_equal(children[["E"]], "H");
    expect_equal(children[["F"]], c("H","I","J"));
    expect_equal(children[["G"]], c("H","E"));
    expect_equal(children[["H"]], character(0));
    expect_equal(children[["I"]], character(0));
    expect_equal(children[["J"]], character(0));
})

test_that("build.ancestors works", {
    g <- make.graph();
    nd <- nodes(g);
    root <- root.node(g);
    lev <- graph.levels(g, root=root);
    anc <- build.ancestors(g);
    anc.tod <- build.ancestors.per.level(g,lev);
    anc.bup <- build.ancestors.bottom.up(g,lev);

    expect_equal(length(anc), 10);
    expect_equal(length(anc.tod), 10);
    expect_equal(length(anc.bup), 10);
    expect_equal(names(anc), nd);
    expect_equal(names(anc.tod), c("A","B","C","G","D","E","F","H","I","J"));
    expect_equal(names(anc.bup), c("J","I","H","F","E","D","G","C","B","A"));
    expect_equal(anc[["A"]], "A");
    expect_equal(anc[["B"]], c("A","B"));
    expect_equal(anc[["C"]], c("A","C"));
    expect_equal(anc[["D"]], c("B","A","D"));
    expect_equal(anc[["E"]], c("A","G","E"));
    expect_equal(anc[["F"]], c("D","B","A","C","F"));
    expect_equal(anc[["G"]], c("A","G"));
    expect_equal(anc[["H"]], c("F","D","B","A","E","G","C","H"));
    expect_equal(anc[["I"]], c("F","D","B","A","C","I"));
    expect_equal(anc[["J"]], c("F","D","B","A","C","J"));
})

test_that("build.descendants works", {
    g <- make.graph();
    nd <- nodes(g);
    root <- root.node(g);
    desc <- build.descendants(g);
    lev <- graph.levels(g, root=root);
    desc.tod <- build.descendants.per.level(g,lev);
    desc.bup <- build.descendants.bottom.up(g,lev);

    expect_equal(length(desc), 10);
    expect_equal(length(desc.tod), 10);
    expect_equal(length(desc.bup), 10);
    expect_equal(names(desc), nd);
    expect_equal(names(desc.tod), c("A","B","C","G","D","E","F","H","I","J"));
    expect_equal(names(desc.bup), c("J","I","H","F","E","D","G","C","B","A"));
    expect_equal(desc[["A"]], c("G","E","H","C","F","J","B","D","I","A"));
    expect_equal(desc[["B"]], c("H","F","J","D","I","B"));
    expect_equal(desc[["C"]], c("H","F","J","I","C"));
    expect_equal(desc[["D"]], c("H","F","J","I","D"));
    expect_equal(desc[["E"]], c("H","E"));
    expect_equal(desc[["F"]], c("H","J","I","F"));
    expect_equal(desc[["G"]], c("E","H","G"));
    expect_equal(desc[["H"]], c("H"));
    expect_equal(desc[["I"]], c("I"));
    expect_equal(desc[["J"]], c("J"));
})

test_that("constraints.matrix works", {
    g <- make.graph();
    child  <- c(2, 3, 7, 5, 4, 6, 6, 8, 6, 9, 8, 8, 9, 10, 8, 5);
    parent <- c(1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 6, 6, 7, 7);
    nodes  <- c("A", "A", "A", "A", "B", "B", "C", "C", "D", "D", "E", "F", "F", "F", "G", "G");
    m <- cbind(child,parent);
    rownames(m) <- nodes;

    expect_equal(constraints.matrix(g), m);
})

test_that("lexicographical.topological.sort works", {
    if (!check.graph()){ skip("graph package cannot be loaded"); }

    g <- make.graph();
    gL1 <- graph::addEdge(from="A",to="A",g);
    gL2 <- graph::addEdge(from="C",to="C",g);
    
    expect_equal(lexicographical.topological.sort(g), c("A","B","C","D","F","G","E","H","I","J"));
    expect_error(lexicographical.topological.sort(gL1), "input graph g is not a dag", fixed=TRUE);
    expect_error(lexicographical.topological.sort(gL2), "input graph g is not a dag", fixed=TRUE);
})

test_that("build.consistent.graph works", {
    if (!check.graph()){ skip("graph package cannot be loaded"); }

    g <- make.graph();
    root <- root.node(g);
    G <- graph::addNode("Z",g);

    expect_s4_class(build.consistent.graph(g, root=root),"graphNEL");
    expect_output(build.consistent.graph(G, root=root), "removed nodes not accessible from root:\\n1 \\t Z");
})

test_that("check.dag.integrity works", { 
    if (!check.graph()){ skip("graph package cannot be loaded"); }

    g <- make.graph();
    root <- root.node(g);
    G <- graph::addNode("Z",g); 
    G <- graph::addEdge(from="Z",to="C",G); 
    G <- graph::addEdge(from="Z",to="Z",G); 

    expect_output(check.dag.integrity(g, root=root), "dag is OK");
    expect_output(check.dag.integrity(G, root=root), 
        "check.GO.integrity: not all nodes accessible from root\nNodes not accessible from root:\nZ", fixed=TRUE);
    expect_error(check.dag.integrity(g, root="R"), "root node not found in g. Insert the root node", fixed=TRUE);
    expect_error(check.dag.integrity(g, root="B"), "the supplied root node is not the right root node of g. Use the function root.node(g) to find the root node of g", fixed=TRUE);
})

test_that("find.leaves works", {
    g <- make.graph();
   
    expect_equal(find.leaves(g), c("H","I","J"));
})

test_that("distances.from.leaves works", {
    g <- make.graph();
    
    dist.leaves <- distances.from.leaves(g);

    nd <- c("A","B","C","D","E","F","G","H","I","J"); 
    dist <- c(2,2,1,1,1,1,1,0,0,0);
    names(dist) <- nd; 

    expect_equal(dist.leaves, dist);
})
