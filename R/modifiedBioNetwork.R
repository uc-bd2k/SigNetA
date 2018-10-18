
#' ModifiedBioNetWork Function
#'
#' This function allows you to analyze genes from a signature using the bionet algorithm from the package "bionet" and uses interactome dataset
#' @param File ? Defaults to NULL.
#' @keywords modifiedbioNetwork
#' @export
#' @examples
#' modifiedBioNetwork() will analyze network for vorinostat signature from ilincs
modifiedBioNetwork<-function(File=NULL,phy=FALSE,layOut=1,package=FALSE,nodeGoData=NULL,edgeGoData=NULL,proteinN=1){

  logic<-File
  ##STRING NETWORK ONE(GET MODULE)##  
  if(proteinN==1)
  {
  #geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
    geninfo<- geneData[which(geneData$GeneID%in%as.character(logic$GeneID)),]
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  pval<-as.numeric(logic$Pvals)
  pval<- -log10(pval)
  names(pval)<-geneLabels
  
  logFC<-as.numeric(logic$coefficients)
  names(logFC)<-geneLabels
  
  
  ##subnet <- subNetwork(geneLabels, interactome)
 ## subnet <- rmSelfLoops(subnet)
  
  ####Modified function used#######
#   load(system.file("extdata", "weightedGraphStringPPI_10.rda", package = "SigNetA"))
#   load(system.file("extdata", "lincscp_1.rda", package = "SigNetA"))
  #ppiGW.copy <- delete.edges(ppiGW, which(E(ppiGW)$weight <=0.7))
  #ppi <- rmSelfLoops(ppiGW.copy)
  ppi<-igraph::simplify(ppiGW,remove.loops = TRUE,remove.multiple = FALSE)
  ppi=decompose.graph(ppi)[[1]]
  
  ###Identify module using FastHeinz algorithm, nsize is fixed to 30 nodes
  
 
 
  names(pval)<-logic$GeneID
  

  module=modules_FastHeinz(subnet=ppi, data_vector=pval)
 
  
  # source("rashidplotmodule.R")
 # pdf("wor.pdf")
  colorNet<-plotmodule2(module, scores =  V(module)$score, diff.expr = logFC)

  module<-igraph.to.graphNEL(colorNet$n) #STRING
  
 # dev.off()

  }
  else if(proteinN==2){
    
    #geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
    geninfo<- geneData[which(geneData$GeneID%in%as.character(logic$GeneID)),]
    geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
    pval<-as.numeric(logic$Pvals)
    pval<- -log10(pval)
    names(pval)<-geneLabels
    
    logFC<-as.numeric(logic$coefficients)
    names(logFC)<-geneLabels
    
    
    
    ####Modified function used#######
    
    ppi<- rmSelfLoops(interactome)
    
    #ppi=decompose.graph(ppi)[[1]] #get the largest subgraph
    
    ###Identify module using FastHeinz algorithm, nsize is fixed to 30 nodes
    
    
    
    names(pval)<-logic$GeneID
    
    
    module=modules_FastHeinz(subnet=ppi, data_vector=pval)
   
    #pdf("wor.pdf")
    colorNet<-plotmodule2(module, scores =  V(module)$score, diff.expr = logFC)
    
    
    module<-igraph.to.graphNEL(colorNet$n) #STRING
  
   # dev.off()
    
    
    
    
    
  }
  
  
  ##END...STRING NETWORK  ONE(GET MODULE##  
  
  ## IGRAPH LAYOUTS 
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
  #geninfo2<-geneInfoFromPortals(geneList=as.character(id),symbol=T,names=F) #STRING
  geninfo2<- geneData[which(geneData$GeneID%in%as.character(id)),]
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
    return (((x*2 - min(x)) / (max(x) - min(x)))*50)
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
  
  ##STRING TWO(ADDING EDGE VALUES)#
  
  if(proteinN=="1"){
    
    
    edgeVisData$title<-"ppi"
    #edgeVisDataMod<-merge(edgeVisData,s,by.x = c("from","to"),by.y =c("a","b"),all.x = TRUE)
    edgeVisDataMod<-igraph::as_data_frame(igraph.from.graphNEL(module),what="edges")
    for(i in 1:nrow(edgeVisDataMod)){
      edgeVisDataMod$title[i]<-paste0("<p><b>Neighborhood score:</b>",edgeVisDataMod$f.neighborhood[i],"</p><b>Fusion score:</b>",edgeVisDataMod$f.fusion[i],"</p><b>Cooccurence score:</b>",edgeVisDataMod$f.cooccurence[i],"</p><b>Coexpression score:</b>",edgeVisDataMod$f.coexpression[i],"</p><b>Experimental score:</b>",edgeVisDataMod$f.experimental[i],"</p><b>Database score:</b>",edgeVisDataMod$f.database[i],"</p><b>Textmining score:</b>",edgeVisDataMod$f.textmining[i],"</p><b>Combined score:</b>",edgeVisDataMod$f.combined_score[i],"</p>")
    }
    edgeVisData<-edgeVisDataMod
  }
  
  ##END...STRING TWO (ADDING EDGE VALUES)##
  if(is.null(nodeGoData) & is.null(edgeGoData))
  {
    
    visObj<- visNetwork(nodeVisData, edgeVisData,height="800px",width="900px")
    
  }
  else{
    visObj<- visNetwork(nodeGoData, edgeGoData,height="800px",width="900px")
    # visObj<-enrichObj
  }


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
  
  return(list("networkObj"=visObj,"networkData"=statNet,"edgeData"=edgeVisData,"nodeData"=nodeVisData))
  
  
} #end of bionet algorithm