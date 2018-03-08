#' D3 Top Hundred Network Function
#'
#' This function allows you to analyze top 100 genes from a signature using the interactome dataset using D3.js
#' @param File ? Defaults to NULL.
#' @keywords topHundredNetwork
#' @export
#' @examples
#' D3layout()  will analyze network for vorinostat signature from ilincs visualized using D3.js
D3layout<-function(File=NULL,upload3=NULL){
  library(igraph)
  library(BioNet)
  library(DLBCL)
  data(interactome)
  library(networkD3)

  if(!is.null(File))
  {
    if(is.null(upload3)){
      File<-File
    }
    else{
      File<-read.csv(file=File,sep='\t')
    }
  }
  else{
    File<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
  }
  sortedFile<-sortNetwork(File)
  logic<-sortedFile
  
  geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  
  #load("/opt/raid10/genomics/mario/ucGithub/netlincs/scratch/data/weightedGraphStringPPI_10.rda")
  #ppiGW.copy <- delete.edges(ppiGW, which(E(ppiGW)$weight <=0.7))
  
  
  subnet <- subNetwork(geneLabels, interactome,neighbors = "none") 
  #subnet <- subNetwork(geneLabels, ppiGW,neighbors = "none") 
  
  subnet <- rmSelfLoops(subnet)
  
  logFC<-as.numeric(logic$coefficients)
  names(logFC)<-geneLabels
  module<-subnet
  colorNet<-plotmodule2(module, diff.expr = logFC)
  
  dev.off();
  
  
  
  
  #library(rcytoscapejs)
  
  
  id <- nodes(module)
  name <- id
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  
  nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
  nodeData$shape <- "ellipse"  #default shape
  nodeData$href <- paste0("http://www.ncbi.nlm.nih.gov/gene/",gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name))))
  nodeData$geneID<-gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name)))
  nodeData$name<-sub(" *\\(.*", "", nodeData$name)
  nodeData$Diff_Exp="none"
  for(i in 1:length(name)){
    nodeData[i,3]<-colorNet$c[i];
    nodeData[i,7]<-colorNet$d[i]
  }
  
  statNet<<-nodeData
  
  ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]])))
  
  
  source<-unlist(lapply(1:length(ltn),function(x) rep(id[x],ltn[x])))
  target<-unlist(lapply(edgeL(module), function(x) id[unlist(x)]))
  networkData<-data.frame(source,target)
  data(MisLinks)
  data(MisNodes)
  pdf("d3.pdf")

 simpleNetwork(networkData)
  dev.off()

  
  
}





