
#' RWR Function
#'
#' This function allows you to analyze genes from a signature using the random walk with restart algorithm and uses STRING dataset
#' @param File ? Defaults to NULL.
#' @keywords RWR
#' @export
#' @examples
#' RWR() will analyze network for vorinostat signature from ilincs
RWR<-function(File=NULL,phy=FALSE,layOut=1,package=FALSE,nodeGoData=NULL,edgeGoData=NULL,proteinN=1){
 
   
  logic<-File
  
  
  
  
##STRING NETWORK ONE(GET MODULE)##  
if(proteinN==1){
  

  geninfo<-data.frame("GeneID"=logic$GeneID,"Symbol"=logic$GeneNames)
 
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
 
  pval<-as.numeric(logic$Pvals)
  pval<- -log10(pval)
  names(pval)<-geneLabels
  
  logFC<-as.numeric(logic$coefficients)
  
  names(logFC)<-geneLabels
  
  

  
 
  
  
  
  
 
  ####Modified function used#######
 
  ppi<-igraph::simplify(ppiGW,remove.loops = TRUE,remove.multiple = FALSE)
  ppi=decompose.graph(ppi)[[1]] #get the largest subgraph
  
  ###Identify module using FastHeinz algorithm, nsize is fixed to 30 nodes
  
  
  
  names(pval)<-logic$GeneID
  
  #added code
  module=modules_RWR_TopScores(subnet=ppi, data_vector=pval,damping_factor=0.8, nseeds=10)

  #CONSTRUCT module only dataframe##
  modIds<-V(module)$name
  geninfo<- geneData[which(as.character(geneData$GeneID) %in% as.character(modIds)),]
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  geninfo$diffExp<-NA
  
  for(i in 1:nrow(geninfo)){
    if(geninfo$GeneID[i] %in% logic$GeneID){
      geninfo$diffExp[i]<-logic$coefficients[match(geninfo$GeneID[i],logic$GeneID)]      
      
    }
    
    }
  logFC<-as.numeric(geninfo$diffExp)
  
  names(logFC)<-geneLabels
  
  
  
  #END of CONSTRUCT module only dataframe##
 

 V(module)$score[is.na(V(module)$score)]<- 0.01




#added code end

  colorNet<-plotmodule2(module, scores =  V(module)$score, diff.expr = logFC)

  
  module<-igraph.to.graphNEL(colorNet$n) #STRING
  

  
  
}
  
else if(proteinN==2){
    

  
  geninfo<-data.frame("GeneID"=logic$GeneID,"Symbol"=logic$GeneNames)
  
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  
  pval<-as.numeric(logic$Pvals)
  pval<- -log10(pval)
  names(pval)<-geneLabels
  
  logFC<-as.numeric(logic$coefficients)
  
  names(logFC)<-geneLabels
  
 
  
  
  
  ####Modified function used#######
  
  ppi<- rmSelfLoops(interactome)

  
  ###Identify module using FastHeinz algorithm, nsize is fixed to 30 nodes
  
  
  
  names(pval)<-logic$GeneID
  
  
  module=modules_RWR_TopScores(subnet=ppi, data_vector=pval, damping_factor=0.8, nseeds=10)
  #CONSTRUCT module only dataframe##
  modIds<-V(module)$name
  geninfo<- geneData[which(as.character(geneData$GeneID) %in% as.character(modIds)),]
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  geninfo$diffExp<-NA
  
  for(i in 1:nrow(geninfo)){
    if(geninfo$GeneID[i] %in% logic$GeneID){
      geninfo$diffExp[i]<-logic$coefficients[match(geninfo$GeneID[i],logic$GeneID)]      
      
    }
    
  }
  logFC<-as.numeric(geninfo$diffExp)
  
  names(logFC)<-geneLabels

  
  
  #END of CONSTRUCT module only dataframe##
  
  #pdf("wor.pdf")
  V(module)$score[is.na(V(module)$score)]<- 0.01
  colorNet<-plotmodule2(module, scores =  V(module)$score, diff.expr = logFC)

  
  module<-igraph.to.graphNEL(colorNet$n) #STRING
  #dev.off()
    
    
    
    
    
  }
##END...STRING NETWORK  ONE(GET MODULE)##  

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
  
  
  ##END..IGRAPH LAYOUT##
  
  
  id <- nodes(module)
  name <- id
  label<-id #visNetwork and interactome
 
 
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
    #return (((x*2 - min(x)) / (max(x) - min(x)))*50)
    return (   20*((x - min(x)) / (max(x) - min(x))) + 30  )
  }
  
  #scoreLength<-max(colorNet$sc)*2 
  
  
  for(i in 1:length(label)){
    nodeVisData[i,3]<-colorNet$c[i];
    nodeVisData[i,6]<-paste0("<p><b>Gene name:</b>",statNet$name[i],"</p><br><p><b>Gene ID:</b>",statNet$geneID[i],"</p><br><p><b>Differential Expression:</b>",statNet$Diff_Exp[i],"</p><br><p><b>Score:</b>",statNet$sc[i],"</p><p><b>NCBI link:</b><a href='",statNet$href[i],"' target='_blank'>",statNet$href[i],"</a></p>")
  
    nodeVisData[i,8]<-colorNet$sc[i]

 
    
  }

  nodeVisData<-data.frame(nodeVisData[1:7], apply(nodeVisData["size"],2, normalize) )
  # for( l in 1:length(nodeVisData$size)){
  #   
  #   if(nodeVisData$size[l]<1)
  #   {
  #     nodeVisData$size[l]<-1
  #   }
  # }
  
  
  
  
  

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
  
  
  
}