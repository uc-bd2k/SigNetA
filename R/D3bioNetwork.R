
#' D3bioNetWork Function
#'
#' This function allows you to analyze genes from a signature using the bionet algorithm from the package "bionet" and uses interactome dataset.Visualization done using D3.js
#' @param File ? Defaults to NULL.
#' @keywords D3bioNetwork
#' @export
#' @examples
#' D3bioNetwork() will analyze network for vorinostat signature from ilincs with D3.js visualization
D3bioNetwork<-function(File=NULL,upload4=NULL){
  library(igraph)
  library(BioNet)
  library(DLBCL)
  data(interactome)
  library(networkD3)
 
  if(!is.null(File))
  {
    if(is.null(upload4)){
      logic<-File
    }
    else{
      logic<-read.csv(file=File,sep='\t')
    }
  }
  else{
    logic<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
  }
  
  
 # source("geneInfoFromPortals.R")
  
  geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  pval<-as.numeric(logic$Pvals)
  
  names(pval)<-geneLabels
  
  logFC<-as.numeric(logic$coefficients)
  names(logFC)<-geneLabels
  
  
  subnet <- subNetwork(geneLabels, interactome)
  subnet <- rmSelfLoops(subnet)
  
  
  system.time( fb <- fitBumModel(pval, plot = FALSE))
  #err2<<-try(scoreNodes(subnet, fb, fdr = 0.1),silent=TRUE)
  #if(class(err2)=="try-error"){
  # output$input_error=renderText("No significant subnetwork generated.Please upload another Signature.")
  # }
  #else{
  #output$input_error=renderText("")
  system.time(scores <- scoreNodes(subnet, fb, fdr = 0.1))
  
  #err<<-try(runFastHeinz(subnet, scores),silent=TRUE)
  # if(class(err) == "try-error"){
  #   
  #   
  #   output$input_error=renderText("No significant subnetwork generated.Please upload another Signature.")
  #   stopifnot(class(err) == "try-error")
  #   
  # }
  
  
  # else{
  #output$input_error=renderText("")
  system.time(module <- runFastHeinz(subnet, scores))
  
  #source("rashidplotmodule.R")
  pdf("wor.pdf")
  colorNet<-plotmodule2(module, scores = scores, diff.expr = logFC)
  dev.off()
  #library(rcytoscapejs)
  
  
  id <- nodes(module)
  name <- id
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  
  nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
  nodeData$shape <- "none"  #default shape
  nodeData$href <- paste0("http://www.ncbi.nlm.nih.gov/gene/",gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name))))
  nodeData$geneID<-gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name)))
  nodeNameEntrez<-nodeData$name
  nodeData$name<-sub(" *\\(.*", "", nodeData$name)
  
  nodeData$Diff_Exp="none"
  nodeData$score="none"
  for(i in 1:length(name)){
    nodeData[i,3]<-colorNet$c[i];
    nodeData[i,7]<-colorNet$d[i]
    nodeData[i,8]<-colorNet$sc[i]
    
  }
  for(i in 1:length(name)){
    if(colorNet$s[i]=="csquare")
      #colorNet$s[i]<-"rectangle"
      colorNet$s[i]<-"ellipse"
    
    else
      colorNet$s[i]<-"ellipse"
    
    nodeData[i,4]<-colorNet$s[i];
    
  }
  statNet<<-nodeData
  
  ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]])))
  source<-unlist(lapply(1:length(ltn),function(x) rep(id[x],ltn[x])))
  target<-unlist(lapply(edgeL(module), function(x) id[unlist(x)]))
  networkData<-data.frame(source,target)
  pdf("d3.pdf")
  
  simpleNetwork(networkData)
  dev.off()
  
  
} #end of bionet algorithm