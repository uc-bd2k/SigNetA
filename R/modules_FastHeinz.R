modules_FastHeinz=function(subnet,data_vector)
  {
    nsize=30
     
    gene_scores=abs(data_vector)[names(data_vector) %in% V(subnet)$name]
    df=data.frame(name=names(gene_scores), gene_scores=gene_scores, rank_gene_scores=length(gene_scores)+1-rank(gene_scores,ties.method ="random"), stringsAsFactors=F)

    positive_nodes=seq(nsize,nsize*3,5)
     
    for (j in 1:length(positive_nodes)) {
      
      idx=positive_nodes[j]+1
      gene_name=df$name[df$rank_gene_scores==idx]   #all ranks are unique
      scores_normalized=gene_scores-df$gene_scores[df$name==gene_name]
  
      scores_other=rep(min(scores_normalized), vcount(subnet)-length(gene_scores))
      names(scores_other)=setdiff(V(subnet)$name, names(gene_scores))
      scores=c(scores_normalized,scores_other)
      module_test <- runFastHeinz(subnet, scores)
      nsize_test <- vcount(module_test)
    if (nsize_test >= nsize) {
      break
      }
    }
    V(module_test)$score=data_vector[V(module_test)$name]
    module_test
  }









