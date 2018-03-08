#' Top Hundred Network Function
#'
#' This function allows you to analyze top 100 genes from a signature using the interactome dataset
#' @param File ? Defaults to NULL.
#' @keywords topHundredNetwork
#' @export
#' @examples
#' topHundredNetwork()  will analyze network for vorinostat signature from ilincs
topHundredNetwork<-function(File=NULL,upload1=NULL,layOut=1,proteinN=2,phy=FALSE,enrich=NULL,package=FALSE){
   
  library(igraph)
  library(BioNet)
  library(DLBCL)
  data(interactome)
  library(networkD3)
  library(visNetwork)

  
  
  if(!is.null(File))
  {  if(is.null(upload1)){
       File<-File
       }
   
       # colnames(File)<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
      
      else{
    File<-read.csv(file=File,sep='\t')
    colnames(File)<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
  
      
    }
  }
  else{
    File<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')

  }
  
  sortedFile<-sortNetwork(File)
  logic<-sortedFile
  print("1")
  if(proteinN=="1"){
  geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F) #interactome
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep="")) #interactome
  }
 print("2")
  if(proteinN=="2"){
  
     
  load(system.file("extdata", "weightedGraphStringPPI_10.rda", package = "SigNetA"))#STRING
    #stringToEntrez<-read.table(file=system.file("extdata", "entrez_gene_id.vs.string.v10.28042015.tsv", package = "SigNetA"),sep='\t')
   
  #  print("ENTREZ TO STRING")
     #print(head(stringToEntrez$V2))
   # a<-stringToEntrez$V2[match(logic$GeneID,stringToEntrez$V1)]
   # print("This is list of a:")
   # print(a)
     ppiGW.copy <- delete.edges(ppiGW, which(E(ppiGW)$weight <=0.7))#STRING
    # subnet <- subNetwork(a, ppiGW.copy,neighbors = "none")
     subnet <- subNetwork(logic$GeneID, ppiGW.copy,neighbors = "none")

    
  }
 print("3")
if(proteinN=="1"){
 subnet <- subNetwork(geneLabels, interactome,neighbors = "none")  #interactome

}
  
  subnet <- rmSelfLoops(subnet)
  if(proteinN=="2"){
  geninfo<-geneInfoFromPortals(geneList=as.character(V(subnet)),symbol=T,names=F) #STRING
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))#STRING
  }
  #STRING genelabels
  
  logFC<-as.numeric(logic$coefficients)
  names(logFC)<-geneLabels
  module<-subnet
  

  colorNet<-plotmodule2(module, diff.expr = logFC)
  
  

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
  
  ###dev.off();
  
  if(proteinN=="2"){
  module<-igraph.to.graphNEL(colorNet$n) #STRING
  }
  
  
  conNodes<-function(x){
    if(ltn[x]>0)
       id[x]
    
  }
  
 
 
  print("4")
  id <- nodes(module) #interactome
  ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]]))) #interactome

  modNodes<-unlist(lapply(1:length(ltn),conNodes ))
 
  id<-modNodes
  if(proteinN=="1"){
  name <- id # interactome and rcytoscapejs2
    label<-id
  }
  if(proteinN=="2")
  {
  geninfo2<-geneInfoFromPortals(geneList=as.character(id),symbol=T,names=F) #STRING
  name<-apply(geninfo2,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))#STRING
  label<-apply(geninfo2,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  }
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE) ##RCYTOSCAPEJS2
  nodeVisData<-data.frame(id,label,stringsAsFactors = FALSE)
  id<-nodes(module)
  
 
  
  nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
  nodeData$shape <- "ellipse"  #default shape
  nodeData$href <- paste0("http://www.ncbi.nlm.nih.gov/gene/",gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name))))
  nodeData$geneID<-gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name)))
  nodeData$name<-sub(" *\\(.*", "", nodeData$name)
  nodeData$Diff_Exp="none"
  nodeData$x="none"   #x and y are the columns required for manual layouts. Also, use "preset" as layout mode in rcytoscapejs2.R
  nodeData$y="none"

  for(i in 1:length(name)){
    nodeData[i,3]<-colorNet$c[i];
    nodeData[i,7]<-colorNet$d[i]
    nodeData[i,8]<-(l[i,1]+1)*100
    nodeData[i,9]<-(l[i,2] +1)*100
  }
 # print(nodeData)
  statNet<<-nodeData
  
 
  
  source<-unlist(lapply(1:length(ltn),function(x) rep(id[x],ltn[x])))

  target<-unlist(lapply(edgeL(module), function(x) id[unlist(x)]))
  
  vect<-c()
  for(i in 1:length(target))  #extracting the value from the key value pair
    vect[i]<-target[[i]]
  
  
  
  
  edgeData <- data.frame(source, target, stringsAsFactors=FALSE)
  
  print("5")
#   network <- createCytoscapeJsNetwork(nodeData, edgeData)
# 
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
    
  #} #to remove doubly linked edges 
 
  #IMPORTANT INFO:
  #resultant  module  visualized  in  Cytoscape.    Signicantly  upregu-
  #lated genes are coloured in red, genes that show signicant downregulation are
  #coloured in green, for the contrast B vs.  T cells.  The score of the nodes is shown
  #by the shape of the nodes, circles indicate a positive score, diamonds a negative
  #score
  
  
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
 
 
 # print(nodeVisData["size"])
  nodeVisData<-data.frame(nodeVisData[1:7], apply(nodeVisData["size"],2, normalize) )
  #print(nodeVisData["size"])
  #nodeVisData<-data.frame(nodeVisData[1:7], abs(as.data.frame(scale(nodeVisData["size"])) ))
  #print(nodeVisData["size"])
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
  ##extract STRING edge scores###
  if(proteinN=="2"){
    load(system.file("extdata", "stringToEntrezDetailed.Rdata", package = "SigNetA"))#STRING
    
    edgeVisData$title<-"ppi"
    edgeVisDataMod<-merge(edgeVisData,s,by.x = c("from","to"),by.y =c("a","b"),all.x = TRUE)
    for(i in 1:nrow(edgeVisDataMod)){
      edgeVisDataMod$title[i]<-paste0("<p><b>Neighborhood score:</b>",edgeVisDataMod$f.neighborhood[i],"</p><b>Fusion score:</b>",edgeVisDataMod$f.fusion[i],"</p><b>Cooccurence score:</b>",edgeVisDataMod$f.cooccurence[i],"</p><b>Coexpression score:</b>",edgeVisDataMod$f.coexpression[i],"</p><b>Experimental score:</b>",edgeVisDataMod$f.experimental[i],"</p><b>Database score:</b>",edgeVisDataMod$f.database[i],"</p><b>Textmining score:</b>",edgeVisDataMod$f.textmining[i],"</p><b>Combined score:</b>",edgeVisDataMod$f.combined_score[i],"</p>")
    }
    edgeVisData<-edgeVisDataMod
  }
  
 
 
  print("6")
  ##GO nodes add start ###
  if(!is.null(enrich)){
   
 
       
      
      res<-c()
      for(i in 1:nrow(enrich)){
        res[i]<-as.character(enrich$idSel[i])
      }
      
    
     
     
      addGO<-genesInEnrichedCategories(res, nodeData$geneID, funcCategories = "GO", species = "Hs")
     
      from<-c()
      to<-c()
      nodeIdData<-c()
      for(i in 1:nrow(addGO)){
        
        for(j in 1:length(addGO[i,])){
          
          
        
          if(addGO[i,j]==TRUE){
           
            subId<-substring(names(addGO[j]),2)
            
            
            nodeIdData<-append(nodeIdData,subId)
       
            
          }
          
        }
       
        genein<-geneInfoFromPortals(geneList=nodeIdData,symbol=T,names=F) #interactome
        geneLab<-apply(genein,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep="")) #interactome
        GOnodeName<-sub(" *\\(.*", "", geneLab)
        for(k in 1:length(GOnodeName)){
          from<-append(from,addGO[[i,1]])
          to<-append(to,GOnodeName[k])
        }
      
       nodeIdData<-NULL
        
      }
      GOdataframe<-data.frame(from,to)
     
      forGo<-nrow(edgeVisData)
      
    
  edgeVisData<-data.frame(from=edgeVisData$from,to=edgeVisData$to)
  
      edgeVisData<-merge(edgeVisData,GOdataframe,all=TRUE)

     
  
 
      
      l<-nrow(nodeVisData)
      
      for(i in 1:nrow(enrich)){
        nodeVisData[l+i,1]<-as.character(enrich$idSel[i])
        nodeVisData[l+i,2]<-as.character(enrich$idSel[i])
        nodeVisData[l+i,3]<-"#ffff00 "
        nodeVisData[l+i,4]<-2
        nodeVisData[l+i,5]<-"black"
        nodeVisData[l+i,6]<-paste0("<p>",enrich$nameSel[i],"</p>")
        nodeVisData[l+i,7]<-"square"
        nodeVisData[l+i,8]<-20
        
      }
      
      
    #}
  }
  
  
  ##GO nodes add stop##
 visObj<- visNetwork(nodeVisData, edgeVisData,height="800px",width="900px")
 
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

 #print("works")
 
 #if(package==TRUE)
 #{
 #  visNetwork(nodeVisData, edgeVisData,height="800px",width="900px")
# }
#print("not work")
}





