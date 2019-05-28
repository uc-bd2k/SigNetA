#' Top Hundred Network Function
#'
#' This function allows you to analyze top 100 genes from a signature using the interactome dataset
#' @param File ? Defaults to NULL.
#' @keywords topHundredNetwork
#' @export
#' @examples
#' topHundredNetwork()  will analyze network for vorinostat signature from ilincs
topHundredNetwork<-function(File=NULL,phy=FALSE,layOut=1,package=FALSE,nodeGoData=NULL,edgeGoData=NULL,proteinN=1){
   

  
  sortedFile<-sortNetwork(File)
  logic<-sortedFile
  
  

  
  
  ##STRING NETWORK ONE(GET MODULE)##  
  if(proteinN==1){
    
   
    
    
    
    ####Modified function used#######
  
    ppi<-igraph::simplify(ppiGW,remove.loops = TRUE,remove.multiple = FALSE)

    
    ###Identify module using FastHeinz algorithm, nsize is fixed to 30 nodes
    
    
    
   # names(pval)<-logic$GeneID
    
    #module=modules_RWR_TopScores(subnet=ppi, data_vector=pval, damping_factor=0.8, nseeds=10)
   
    module <- subNetwork(as.character(logic$GeneID), ppi,neighbors = "none")
   
    
    #CONSTRUCT module only dataframe##

   modIds<-V(module)$name

   geninfo<- geneData[which(as.character(geneData$GeneID) %in% as.character(modIds)),]
    
    geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
    geninfo$diffExp<-NA
    geninfo$pval<-NA
  
    
    for(i in 1:nrow(geninfo)){
     
      if(geninfo$GeneID[i] %in% logic$GeneID){
        geninfo$diffExp[i]<-logic$coefficients[match(geninfo$GeneID[i],logic$GeneID)]     
        geninfo$pval[i]<-logic$Pvals[match(geninfo$GeneID[i],logic$GeneID)]  
      
      }
      
    }
    logFC<-as.numeric(geninfo$diffExp)
    
    names(logFC)<-geneLabels
   
  
    pval<-as.numeric(geninfo$pval)
  
    pval<- -log10(pval)
    names(pval)<-geneLabels
  
    
    
    #END of CONSTRUCT module only dataframe##
    
   
    
    colorNet<-plotmodule2(module, scores =  pval, diff.expr = logFC)

    
   
    
    colorNetframe<-data.frame(V(colorNet$n)$name,colorNet$sc,colorNet$d)
    module<-igraph.to.graphNEL(colorNet$n) #STRING
    
   
    
   # dev.off()
    
    
  }
  
  else if(proteinN==2){
    
  
    
   # ppi<- rmSelfLoops(interactome)
    ppi<-interactome
   
    
    #ppi=decompose.graph(ppi)[[1]] #get the largest subgraph
    
    ###Identify module using FastHeinz algorithm, nsize is fixed to 30 nodes
    
    
    
   # names(pval)<-logic$GeneID
    
    

    
    module <- subNetwork(as.character(logic$GeneID), ppi,neighbors = "none")
 
    
    modIds<-V(module)$name
    
    
    geninfo<- geneData[which(as.character(geneData$GeneID) %in% as.character(modIds)),]

    geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
    geninfo$diffExp<-NA
    geninfo$pval<-NA
   
    
    for(i in 1:nrow(geninfo)){
      
      if(geninfo$GeneID[i] %in% logic$GeneID){
        geninfo$diffExp[i]<-logic$coefficients[match(geninfo$GeneID[i],logic$GeneID)]     
        geninfo$pval[i]<-logic$Pvals[match(geninfo$GeneID[i],logic$GeneID)]  
        
      }
      
    }
    logFC<-as.numeric(geninfo$diffExp)
    
    names(logFC)<-geneLabels
   
   
    pval<-as.numeric(geninfo$pval)
    
    pval<- -log10(pval)
    names(pval)<-geneLabels
   
    

    colorNet<-plotmodule2(module, scores =  pval, diff.expr = logFC)
    
    
    module<-igraph.to.graphNEL(colorNet$n) #STRING
    
   # dev.off()
    
    
    
    
    
  }
  ##END...STRING NETWORK  ONE(GET MODULE)##  
  
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
  
  ###dev.off();
  
  
  
  
  conNodes<-function(x){
    
    if(ltn[x]>0)
      id[x]
    
  }
  
  
  
 
  
  id <- nodes(module) #interactome

  
  ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]]))) #interactome
  
  modNodes<-unlist(lapply(1:length(ltn),conNodes )) # only get connected nodes and no single nodes
  
  id<-modNodes
  
  #TEST#
  id <- nodes(module)
  name <- id
  label<-id #visNetwork and interactome
  #END TEST#
  
  

  
 
  geninfo2<- geneData[which(geneData$GeneID%in%as.character(id)),]
  name<-apply(geninfo2,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))#STRING
  label<-apply(geninfo2,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))

  nodeData <- data.frame(id, name, stringsAsFactors=FALSE) ##RCYTOSCAPEJS2
  nodeVisData<-data.frame(id,label,stringsAsFactors = FALSE)
  id<-nodes(module)
  

  for(i in 1:nrow(geninfo2)){
    
    if(geninfo2$GeneID[i] %in% logic$GeneID){
      colorNet$d[i]<-logic$coefficients[match(geninfo2$GeneID[i],logic$GeneID)]     
      colorNet$sc[i]<-(-log10(logic$Pvals[match(geninfo2$GeneID[i],logic$GeneID)]  ))
      
    }
    
  }
  
  
  nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
  nodeData$shape <- "ellipse"  #default shape
  nodeData$href <- paste0("http://www.ncbi.nlm.nih.gov/gene/",gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name))))
  nodeData$geneID<-gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name)))
  nodeData$name<-sub(" *\\(.*", "", nodeData$name)
  nodeData$Diff_Exp="none"
  nodeData$score="none"
  nodeData$x="none"   #x and y are the columns required for manual layouts. Also, use "preset" as layout mode in rcytoscapejs2.R
  nodeData$y="none"
  
  for(i in 1:length(name)){
    nodeData[i,3]<-colorNet$c[i];
    nodeData[i,7]<-colorNet$d[i]
    nodeData[i,8]<-colorNet$sc[i]
    nodeData[i,9]<-(l[i,1]+1)*100
    nodeData[i,10]<-(l[i,2] +1)*100
  }
  
 
  statNet<<-nodeData

  #statNet$df_data<<-nodeData
  
  
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
  
  
  
  for(i in 1:length(label)){
    nodeVisData[i,3]<-colorNet$c[i];
    nodeVisData[i,6]<-paste0("<p><b>Gene name:</b>",statNet$name[i],"</p><br><p><b>Gene ID:</b>",statNet$geneID[i],"</p><br><p><b>Differential Expression:</b>",statNet$Diff_Exp[i],"</p><br><p><b>Score:</b>",statNet$score[i],"</p><p><b>NCBI link:</b><a href='",statNet$href[i],"' target='_blank'>",statNet$href[i],"</a></p>")
    # nodeVisData[i,8]<-(scoreLength - colorNet$sc[i])
    nodeVisData[i,8]<-colorNet$sc[i]
    
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
    # edgeVisDataMod<-merge(edgeVisData,s,by.x = c("from","to"),by.y =c("a","b"),all.x = TRUE)
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
    #print("physics active")
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





