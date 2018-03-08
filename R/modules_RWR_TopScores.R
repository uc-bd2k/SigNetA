modules_RWR_TopScores=function(subnet,data_vector,damping_factor,nseeds)
{
  library(igraph)
  library(BioNet)
  nsize=30
  
  input_scores=abs(data_vector)[names(data_vector) %in% V(subnet)$name]
  topgenes=names(input_scores[order(input_scores, decreasing=T)][1:nseeds])
  
  input_scores[names(input_scores)]=0
  input_scores[topgenes]=1
  
  scores_other=rep(min(input_scores), vcount(subnet)-length(input_scores))
  names(scores_other)=setdiff(V(subnet)$name, names(input_scores))
  
  p0=c(input_scores,scores_other)
  p0_sum2one=p0/sum(p0)
  
  junk3=  page_rank(subnet, algo = "prpack", vids = V(subnet), directed = TRUE, damping = damping_factor, 
                    personalized = p0_sum2one[V(subnet)$name], weights = NA, options = NULL)
  gene_scores=junk3[[1]]
  
  positive_nodes=seq(nsize,nsize*10,5)
  
  for (j in 1:length(positive_nodes)) {
    #      print(j)
    idx=positive_nodes[j]
    gene_scores_ordered=gene_scores[order(gene_scores, decreasing=T)]
    geneLabels=names(gene_scores_ordered[1:idx])
    
    modules <- subNetwork(geneLabels, subnet)
    modules =decompose.graph(modules)
    largest <- which.max(sapply(modules, vcount))
    
    module_test=modules[[largest]]
    nsize_test <- vcount(module_test)
    print(nsize_test)
    if (nsize_test >= nsize) {
      break
    }
  }
  V(module_test)$score=data_vector[V(module_test)$name]
  module_test
}