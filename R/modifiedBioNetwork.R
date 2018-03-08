
#' ModifiedBioNetWork Function
#'
#' This function allows you to analyze genes from a signature using the bionet algorithm from the package "bionet" and uses interactome dataset
#' @param File ? Defaults to NULL.
#' @keywords modifiedbioNetwork
#' @export
#' @examples
#' modifiedBioNetwork() will analyze network for vorinostat signature from ilincs
modifiedBioNetwork<-function(File=NULL,upload4=NULL,phy=FALSE,layOut=1,package=FALSE){
  library(igraph)
  library(BioNet)
  library(DLBCL)
  data(interactome)
  library(visNetwork)
  
  if(!is.null(File))
  {  
    if(is.null(upload4)){
      logic<-File
    }
    else{
      logic<-read.csv(file=File,sep='\t')
      colnames(logic)<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
    }
  }
  else{
    logic<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
  }
  
  
  #source("geneInfoFromPortals.R")
  
  geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  pval<-as.numeric(logic$Pvals)
  pval<- -log10(pval)
  names(pval)<-geneLabels
  
  logFC<-as.numeric(logic$coefficients)
  names(logFC)<-geneLabels
  
  
  ##subnet <- subNetwork(geneLabels, interactome)
 ## subnet <- rmSelfLoops(subnet)
  
  ####Modified function used#######
  load(system.file("extdata", "weightedGraphStringPPI_10.rda", package = "SigNetA"))
  load(system.file("extdata", "lincscp_1.rda", package = "SigNetA"))
  ppiGW.copy <- delete.edges(ppiGW, which(E(ppiGW)$weight <=0.7))
  ppi <- rmSelfLoops(ppiGW.copy)
  ppi=decompose.graph(ppi)[[1]]
  
  ###Identify module using FastHeinz algorithm, nsize is fixed to 30 nodes
  
 
 
  names(pval)<-logic$GeneID
  
  #module=modules_FastHeinz(subnet=ppi, data_vector=lincscp_1)
  module=modules_FastHeinz(subnet=ppi, data_vector=pval)
  #module=
  ##modified ends##
  
  ##system.time( fb <- fitBumModel(pval, plot = FALSE))
  #err2<<-try(scoreNodes(subnet, fb, fdr = 0.1),silent=TRUE)
  #if(class(err2)=="try-error"){
  # output$input_error=renderText("No significant subnetwork generated.Please upload another Signature.")
  # }
  #else{
  #output$input_error=renderText("")
  ##system.time(scores <- scoreNodes(subnet, fb, fdr = 0.1))
  
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
  ##system.time(module <- runFastHeinz(subnet, scores))
  
  # source("rashidplotmodule.R")
  pdf("wor.pdf")
  colorNet<-plotmodule2(module, scores =  V(module)$score, diff.expr = logFC)

  module<-igraph.to.graphNEL(colorNet$n) #STRING
  dev.off()
  #library(rcytoscapejs)
  ## IGRAPH LAYOUTS FOR RCYTOSCAPEJS2
  if(layOut=="1"){
    l<-layout_with_fr(colorNet$n)
    visLay<-"layout_with_fr"
  }
  
  else if(layOut=="3"){
    l<-layout_on_grid(colorNet$n)
    visLay<-"layout_on_grid"
    
  }
  else if(layOut=="4"){
    l<-layout_with_kk(colorNet$n)
    visLay<-"layout_with_kk"
  }
  else if(layOut=="5"){
    l<-layout_on_sphere(colorNet$n)
    visLay<-"layout_on_sphere"
  }
  else if(layOut=="6"){
    l<-layout_with_graphopt(colorNet$n)
    visLay<-"layout_with_graphopt"
  }
  
 
  id <- nodes(module)
  name <- id
  label<-id #visNetwork and interactome
  geninfo2<-geneInfoFromPortals(geneList=as.character(id),symbol=T,names=F) #STRING
  name<-apply(geninfo2,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))#STRING
  label<-apply(geninfo2,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))#STRING
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  nodeVisData<-data.frame(id,label,stringsAsFactors = FALSE) #visNetwork
  
 
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
    nodeData[i,7]<-colorNet$d[i];
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
  vect<-c()
  for(i in 1:length(target))  #extracting the value from the key value pair
    vect[i]<-target[[i]]
  
  
  
  
  edgeData <- data.frame(source, target, stringsAsFactors=FALSE)
 
  
#   network <- createCytoscapeJsNetwork(nodeData, edgeData)
#   for(i in 1:length(target)){
#     
#     network$edges[[i]]$data$edgeTargetShape="none"  #making undirected graphss
#     
#   }
#   for(i in 1:length(target)){
#     for(j in i:length(target)){
#       if(network$edges[[i]]$data$source == network$edges[[j]]$data$target)
#         network$edges[[j]]$data$target= "none"
#       
#     }
#     
#   }
  ##VisNetwork###########
  
  nodeVisData$color.background<-rep("blue",nrow(nodeData))
  
  nodeVisData$borderWidth <- 2
  nodeVisData$color.border <- "black"
  nodeVisData$title<-"<p>Hello world</p>"
  nodeVisData$shape<-"dot"
  nodeVisData$size<-0
  nodeVisData$id<-sub(" *\\(.*", "", nodeVisData$id)
  nodeVisData$label<-sub(" *\\(.*", "", nodeVisData$label)
  
  normalize <- function(x) {
    return (((x - min(x)) / (max(x) - min(x)))*50)
  }
  

  for(i in 1:length(label)){
    nodeVisData[i,3]<-colorNet$c[i];
    nodeVisData[i,6]<-paste0("<p><b>Gene name:</b>",statNet$name[i],"</p><br><p><b>Gene ID:</b>",statNet$geneID[i],"</p><br><p><b>Differential Expression:</b>",statNet$Diff_Exp[i],"</p><p><b>NCBI link:</b><a href='",statNet$href[i],"' target='_blank'>",statNet$href[i],"</a></p>")
    if(colorNet$d[i]<0)
    {
      nodeVisData[i,8]<-colorNet$d[i] * -1 
      
    }
    else{
      
      nodeVisData[i,8]<-colorNet$d[i]
    }
    
  }
  

  nodeVisData<-data.frame(nodeVisData[1:7], apply(nodeVisData["size"],2, normalize) )
  for( l in 1:length(nodeVisData$size)){
    
    if(nodeVisData$size[l]<1)
    {
      nodeVisData$size[l]<-1
    }
  }
  
 
  
  
  ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]])))
  
  sourceVis<-unlist(lapply(1:length(ltn),function(x) rep(id[x],ltn[x])))
  targetVis<-unlist(lapply(edgeL(module), function(x) id[unlist(x)]))
  vect<-c()
  for(i in 1:length(targetVis))  #extracting the value from the key value pair
    vect[i]<-targetVis[[i]]
  
  sourceVis<-sub(" *\\(.*", "", sourceVis)
  targetVis<-sub(" *\\(.*", "", targetVis)
  
  
  edgeVisData <- data.frame(from=sourceVis, to=targetVis, stringsAsFactors=FALSE)

  for (i in 1:nrow(edgeVisData))
  {
    edgeVisData[i, ] = sort(edgeVisData[i, ])
  }
  edgeVisData = edgeVisData[!duplicated(edgeVisData),]
  
  visObj<- visNetwork(nodeVisData, edgeVisData,height="800px",width="900px")

  visObj<-visExport(visObj,type ="png", name="network",float="right")
  
  if(phy){
   # print("physics active")
  }
  else{
    visObj<-visIgraphLayout(visObj,layout=visLay)
  }
  
  if(package==TRUE)
  {
    visNetwork(nodeVisData, edgeVisData,height="800px",width="900px")
  }
  else{
    visObj<-visInteraction(visObj,navigationButtons = TRUE)
    
    
    visObj<-visOptions(visObj,manipulation = TRUE)
  }
  
  #plotInput(network,network$nodes,network$edges)
  #rcytoscapejs2(network$nodes, network$edges,width=1500,height=1800, layout="spread", showPanzoom=TRUE)
  #}
  
  ###END FILE GENERATE UPLOAD NETWORK###
  #}
  
  
} #end of bionet algorithm