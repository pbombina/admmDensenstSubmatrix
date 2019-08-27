install.packages("GGally")
library(GGally)

library(network)
library(sna)
library(ggplot2)

ggnet2(new)



g=graph.adjacency(new,mode="undirected",weighted=NULL,diag=FALSE)
plot.igraph(g)

#Use this
quartz()
gplot(new,gmode="graph")


#OR
library(igraph)


#Choose your favorite algorithm to find communities.  The algorithm below is great for large networks but only works with undirected graphs
c_g <- fastgreedy.community(g)

#Collapse the graph by communities.  This insight is due to this post http://stackoverflow.com/questions/35000554/collapsing-graph-by-clusters-in-igraph/35000823#35000823

res_g <- simplify(contract(g, membership(c_g)))
plot(g, margin = -.5)


#OR even better
plot(simplify(g), vertex.size= 0.01,edge.arrow.size=0.001,vertex.label.cex = 0.75,vertex.label.color = "black"  ,vertex.frame.color = adjustcolor("white", alpha.f = 0),vertex.color = adjustcolor("white", alpha.f = 0),edge.color=adjustcolor(1, alpha.f = 0.15),display.isolates=FALSE,vertex.label=ifelse(page_rank(g)$vector > 0.1 , "important nodes", NA))
