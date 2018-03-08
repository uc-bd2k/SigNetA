#' dmGWAS Function
#'
#' This function allows you to analyze genes  using  algorithm mentioned in Ideker et.al. 2002
#' @param File ? Defaults to NULL.
#' @keywords dmGWAS
#' @export
#' @examples
#' dmGWAS()  will analyze network for genes in signature from ilincs using ideker et. al. 2002 algorithm
dmGWAS<-function(File=NULL,upload3=NULL,layOut=1,proteinN=2,phy=FALSE,enrich=NULL,package=FALSE){
library(dmGWAS)
library(DLBCL)
library(visNetwork)
data(interactome)
if(!is.null(File))
{  if(is.null(upload3)){
  File<-File
}
  else{
    File<-read.csv(file=File,sep='\t')
    colnames(logic)<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
  }
}
else{
  File<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')

}
#fileSNP<-read.csv(file="ExtraFiles/dmGWAS/diabetesSNP.csv",sep=",")
#gene.map=SNP2Gene.match(assoc.file="ExtraFiles/dmGWAS/diabetesSNP.csv",snp2gene.file="ExtraFiles/dmGWAS/GenomeWideSNP_5.na30.annot.AffyID", boundary=20)

##dataFile<- read.csv(file="ExtraFiles/dmGWAS/vorinostat.xls",sep="\t")
#GeneNames and Pvals

gene2weight<-data.frame(gene=File$GeneNames,weight=File$Pvals)
sorted<-gene2weight[order(gene2weight$weight),] 

topgenes<-head(sorted,201)

gene2weight<-topgenes

print("stage 1 passed")

i_interactome<-igraph.from.graphNEL(interactome)

ppinetwork<-as.data.frame(get.edgelist(i_interactome))
ppinetwork$V1<-sub(" *\\(.*", "", ppinetwork$V1)
ppinetwork$V2<-sub(" *\\(.*", "", ppinetwork$V2)
print(head(ppinetwork))
#p_val<-c()
#for(i in 1:length(gene2weight$weight))
# { 

#if(gene2weight$weight[i]>=1)
#{   
# gene2weight$weight[i]<-0.999
#append(p_val,gene2weight$weight[i])
#}
# else if(gene2weight$weight[i]<=0){
# gene2weight$weight[i]<-0.0000001
# }
#}
print("stage 2 passed")
res.list = dms(ppinetwork, gene2weight, d=2, r=0.1) #works for 300 genes and d=1

##error while running res.list with different parameters
##100 genes error is Error in dm.result[[k]] : subscript out of bounds
##200 genes error is -Error in identical.idx[[k]] : subscript out of bounds
##201 and above works
##
print("stage 3 passed")
selected = simpleChoose(res.list, top=100, plot=T) 
print("stage 4 passed")

#source("ExtraFiles/dmGWAS/plotmodule2.R")
logFC<-as.numeric(File$coefficients)
names(logFC)<-File$GeneNames
colorNet<-plotmodule2(selected$subnetwork, diff.expr =logFC)
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
module<-igraph.to.graphNEL(selected$subnetwork)

conNodes<-function(x){
  if(ltn[x]>0)
    id[x]
  
}


id <- nodes(module) #interactome
ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]]))) #interactome

modNodes<-unlist(lapply(1:length(ltn),conNodes ))

id<-modNodes

name <- id # interactome and rcytoscapejs2
label<-id

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
  nodeData[i,7]<-colorNet$d[[i]];
  # nodeData[i,8]<-(l[i,1]+1)*100
  # nodeData[i,9]<-(l[i,2] +1)*100
}

statNet<<-nodeData



nodeVisData$color.background<-rep("blue",nrow(nodeVisData))

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


print(label)
for(i in 1:length(label)){
  nodeVisData[i,3]<-colorNet$c[i];
  nodeVisData[i,6]<-paste0("<p><b>Gene name:</b>",statNet$name[i],"</p><br><p><b>Differential Expression:</b>",statNet$Diff_Exp[i])
  # nodeVisData[i,6]<-paste0("GENE NAME");
  if(colorNet$d[[i]]<0)
  {
    nodeVisData[i,8]<-colorNet$d[[i]] * -1 
    
  }
  else{
    
    nodeVisData[i,8]<-colorNet$d[[i]]
  }
  
}

print(nodeVisData)

nodeVisData<-data.frame(nodeVisData[1:7], apply(nodeVisData["size"],2, normalize) )
for( l in 1:length(nodeVisData$size)){
  print(nodeVisData$size[l])
  if(nodeVisData$size[l]<1)
  {
    nodeVisData$size[l]<-1
  }
}

print(nodeVisData)


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
print(visObj)
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


}




