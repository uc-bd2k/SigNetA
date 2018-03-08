
library(DT)
library(shiny)
library(shinyjs)
library(DLBCL)
library(BioNet)
library(igraph)
library(CLEAN)
library(CLEAN.Hs)
library(networkD3)
library(SigNetA)
library(visNetwork)
library(dmGWAS)
#library(RCurl)

data(interactome)
statNet<<-NULL
visImg<<-NULL
exFile<<-NULL
en<<-NULL
tableEnrichment<<-NULL
shinyServer(function(input,output,session){
  observe({
    query<-parseQueryString(session$clientData$url_search)
    
    if (!is.null(query[["File"]])) {
      File<<-query[['File']]
      ##print(File)
      #print("woohoo")
      # File<<- read.delim(file=File)
      # File<<- read.csv(file=File,sep="\t")
      # #print(class(File))
      # #print(data[[GeneNames]])
      path<-paste("http://dev.ilincs.org/tmp/",query[["File"]],sep="") #get file from url look it up
     # path<-paste("/Users/Rashid/Desktop/Rashid/Career/PhD/Research/Events/BD2KAllHandsMeeting/signatures/",query[["File"]],sep="")
      File<<-path
      # File<<- read.csv(file=path,sep="\t")
  #File<<-File[complete.cases(File),]
  
  shinyjs::enable("downloadNetworkImage")
  shinyjs::show("saveImage")
  shinyjs::enable("downloadEnrich")
  shinyjs::show("algorithm")
  shinyjs::show("layout")
  shinyjs::show("PPI")
  shinyjs::show("goAn")
  toggle(condition = File  , selector = "#tabs li a[data-value=cytonet]")
  toggle(condition = File, selector = "#tabs li a[data-value=datTable]")
  toggle(condition = File, selector = "#tabs li a[data-value=irich]")
  
Sys.sleep(2)
  
  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content", anim = TRUE, animType = "fade")
      ##print(path)
      #File<<- read.delim(file=path)
    }
    
    else{
      File<<-""
    }
  })
 observeEvent(input$load_example,{
    
    exFile<<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
  
   shinyjs::enable("downloadNetworkImage")
   shinyjs::show("saveImage")
   shinyjs::enable("downloadEnrich")
   shinyjs::show("algorithm")
   shinyjs::show("layout")
   shinyjs::show("PPI")
   shinyjs::show("goAn")
   toggle(condition = TRUE  , selector = "#tabs li a[data-value=cytonet]")
   toggle(condition = TRUE , selector = "#tabs li a[data-value=datTable]")
   toggle(condition = TRUE , selector = "#tabs li a[data-value=irich]")
   
  })
 #shinyjs::toggle(selector = "span.logo",TRUE,condition="[data-collapsed=true]" )

  #shinyjs::onclick("[role=button]",shinyjs::hide(selector = "[role=button]"))
 
 
  observe({
  #shinyjs::onclick(".skin-yellow .main-header .navbar .sidebar-toggle ",shinyjs::hide(".skin-yellow .main-header .logo "))
    if (is.null(input$file1) || input$file1 == ""  ) {
       print(grepl("xls",File))
      
        if(!(grepl("xls",File))){
          print("grepl")
      shinyjs::disable("downloadNetworkImage")
      shinyjs::hide("saveImage")
      shinyjs::disable("downloadEnrich")
      shinyjs::hide("algorithm")
      shinyjs::hide("layout")
      shinyjs::hide("PPI")
      shinyjs::hide("goAn")
        }
      else{
        
      }
      
    } else {
    print("YES")
      shinyjs::enable("downloadNetworkImage")
      shinyjs::show("saveImage")
      shinyjs::enable("downloadEnrich")
      shinyjs::show("algorithm")
      shinyjs::show("layout")
      shinyjs::show("PPI")
      shinyjs::show("goAn")
      
      
      Sys.sleep(2)
      
      # Hide the loading message when the rest of the server function has executed
      hide(id = "loading-content", anim = TRUE, animType = "fade")
      
    }
    
    if(is.null(File) || File==""){
    toggle(condition = input$file1  , selector = "#tabs li a[data-value=cytonet]")
    toggle(condition = input$file1, selector = "#tabs li a[data-value=datTable]")
    toggle(condition = input$file1, selector = "#tabs li a[data-value=irich]")
    }
    
    
    
    
  })
  observeEvent(input$saveImage, {
    
    session$sendCustomMessage(type="saveImage", message="NULL")
    
  })
 
observe({
  if(input$algorithm==3 || input$algorithm==4){
    showModal(modalDialog(
      title = "Important message",
      "Please wait.This algorithm takes time to generate the subnetwork."
    ))
  }
})
  
  output$clickedNode = renderPrint({
    input$clickedNode
  })
  
  output$connectedNodes = renderPrint({
    input$connectedNodes
  })
  
  ##gdgsgds
  output$enrichmentAnalysis<-DT::renderDataTable({
    if (is.null(input$file1) || input$file1 == "" ) {
      if(!is.null(exFile)){
        
        enrichmentChange()
        
      }
      else{
        paste("No genes analyzed")
      }
        
      
    } else if(!is.null(input$file) || !is.null(input$PPI)){
      
      
     # if(class(err) == "try-error"){ 
        
       # output$input_error2=renderText("No significant enrichment analysis data generated.Please upload another signature.")
      #  stopifnot(class(err) == "try-error")
      #}
     # else{
       # output$input_error2=renderText("")
     
        enrichmentChange()
        
       #this is where enrichment functionality was
        
        
      }
      
      
  #  }
  }
  
  
  )
  
  #enrichmentChange<-reactive({
   enrichmentChange<-function(){
   # if((input$algorithm=="1" && input$layout=="1") || (input$algorithm=="1" && input$layout=="2"))
    if(input$algorithm=="1")
    {
      if(is.null(input$file1)){
        logic<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
      }
     
      else{
      logic<-read.csv(input$file1$datapath,sep="\t")
      }
      
      
      logic<-sortNetwork(logic)
     
      
      geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
      geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
     
      if(input$PPI=="2"){
        subnetGenes<-statNet$id
      }
      else{
      subnetGenes<-gsub("[\\(\\)]", "", regmatches( statNet$id, gregexpr("\\(.*?\\)",  statNet$id)))
      }
      allGenes<-gsub("[\\(\\)]", "", regmatches(geneLabels, gregexpr("\\(.*?\\)", geneLabels)))
   

   
   
      result<-geneListEnrichment(subnetGenes, allGenes, functionalCategories = "GO", species = "Hs", minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, verbose=TRUE)
    
      
      
      
      options(scipen = 999)
      ID<-result$categories
      Name<-result$Description
      PValue<-format(result$FisherPValue,digits=2,scientific = TRUE)
      FDR<-format(result$FisherFDR,digits=2,scientific = TRUE)
      NetworkGenes<-result$nGenesInCategory
      AllGenes<-result$nAllGenesInCategory
      logOR<-round(result$logOR,2)
      tableEnrichment<<-data.frame(ID,Name,PValue,FDR,NetworkGenes,AllGenes,logOR)
     datatable(tableEnrichment, options=list(columnDefs=list(list(targets=4)),pageLength = 7))
    
      
   }
    
   # else if((input$algorithm=="2" && input$layout=="2")||(input$algorithm=="2" && input$layout=="1")){
    else if(input$algorithm=="2" || input$algorithm=="4"){
      if(is.null(input$file1)){
        logic<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
      }
      else{
        logic<-read.csv(input$file1$datapath,sep="\t")
      }
      
      
      
   
      
      geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
      geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
      
      subnetGenes<-gsub("[\\(\\)]", "", regmatches( statNet$id, gregexpr("\\(.*?\\)",  statNet$id)))
      allGenes<-gsub("[\\(\\)]", "", regmatches(geneLabels, gregexpr("\\(.*?\\)", geneLabels)))
     
      result<-geneListEnrichment(subnetGenes, allGenes, functionalCategories = "GO", species = "Hs", minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, verbose=TRUE)
      
   
      
      
      options(scipen = 999)
      ID<-result$categories
      Name<-result$Description
      PValue<-format(result$FisherPValue,digits=2,scientific = TRUE)
      FDR<-format(result$FisherFDR,digits=2,scientific = TRUE)
      NetworkGenes<-result$nGenesInCategory
      AllGenes<-result$nAllGenesInCategory
      logOR<-round(result$logOR,2)
      tableEnrichment<<-data.frame(ID,Name,PValue,FDR,NetworkGenes,AllGenes,logOR)
     
      
    }
    
   } 
  #})
  
 # tableChange<-reactive({
  tableChange<-function(){ 

   # if(input$algorithm=="1" ){
    
      statNet<-statNet[-c(1,3,4)]
      
      colnames(statNet)<-c("geneName","NCBI Information","geneID","Diff_Exp")
      statNet<-statNet[c("geneName","geneID","Diff_Exp","NCBI Information")]
      
      
     
      statNet
   # }
  }
 # })
  
  output$contents <- DT::renderDataTable({
    
    if (is.null(input$file1) || input$file1 == "" ) {
      if(!is.null(exFile) || is.null(input$PPI)){
        
        tableChange()
        
      }
      else{
        paste("No genes analyzed")
      }
    }
    
    else if(!is.null(input$file) || !is.null(input$PPI)){
      
      
      
      tableChange()
      
    }
    
    
  })
 
  
  observeEvent(input$goAn,{
 
    
    idSel<<-as.character(tableEnrichment[input$enrichmentAnalysis_rows_selected,1])
    nameSel<<-as.character(tableEnrichment[input$enrichmentAnalysis_rows_selected,2])
    en<<-data.frame(idSel,nameSel)
    
  
  if(nrow(en)>0){
    
   }
   else{
      
      en<<-NULL
     
    }
    
    
  
    output$plot<-renderVisNetwork({
    
      ret<-topHundredNetwork(input$file1$datapath,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=en,package=FALSE)
      
    })
   
  
  
  })
  output$genenodes <- renderText({ 
    
    if (is.null(input$file1) || input$file1 == "") {
      paste("No genes analyzed")
    } else {
      
      
      newlist<-c()
      
      
      for (i in 1:length(statNet$id)){
        append(statNet$id,"\n",i)
        newlist[i]<-statNet$id[i]
      }
      
      newlist
      
      
      
      
    }
    
    
    
    
  })
  
  
 
   
    
    

  
  output$general_ui <- renderUI({
   if(input$algorithm=="1" && input$layout=="1"){
     #rcytoscapejsOutput("plot", height="800px",width="900px")
     visNetworkOutput("plot")
   }
    else if(input$algorithm=="2" &&  input$layout=="1"){
      #rcytoscapejsOutput("plot", height="800px",width="900px")
      visNetworkOutput("plot")
    }
    
    
    #LAYOUT
     else if(input$algorithm=="1" && input$layout=="2"){
      simpleNetworkOutput("simple")
     }
    else if(input$algorithm=="2" && input$layout=="2"){
      simpleNetworkOutput("simple")
    }
    
    else{
      #rcytoscapejsOutput("plot", height="800px",width="900px")
      visNetworkOutput("plot")
    }
    
    
    
    
  })
  
  enrichAdd<-function(){
    
    if(is.null(en)){
       retu<-NULL;
    }
    else{
      retu<-en
    }
  }
    
  
  datasetInput <- reactive({
    if(input$algorithm=="1" && (!grepl("xls",File))){
      #if(!is.null(File))
      #{ print(File)
     #   ret<-topHundredNetwork(File,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=enrichAdd())
      #}
     # else{
        print(File)
      ret<-topHundredNetwork(input$file1$datapath,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=enrichAdd())
      #}
      
    }
    else if(input$algorithm=="1" && (grepl("xls",File))){
      ret<-topHundredNetwork(File,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=enrichAdd())
    }
   
    
    
    else if(input$algorithm=="2"){
      
      
      ret<-bioNetwork(input$file1$datapath,upload2="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
    }
    
    else if(input$algorithm=="3"){
      ret<-dmGWAS(input$file1$datapath,upload3="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
    }
    
    else if(input$algorithm=="4"){
      modifiedBioNetwork(input$file1$datapath,upload4="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
    }
    else if(input$algorithm=="5"){
      RWR(input$file1$datapath,upload5="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
    }
    
   
    
    
    
    
    
  })
  
  layoutInput<-reactive({
    
    #if(input$layout=="2"){
      if(input$algorithm=="1"){
        #source("D3layout.R")
        ret<-D3layout(input$file1$datapath,upload3="yes")
      }
      else if(input$algorithm=="2"){
      #source("D3bioNetwork.R")
      ret<-D3bioNetwork(input$file1$datapath,upload4="yes")
      }
    }
    
    
 # }
    
  )
  

    
    output$plot<-renderVisNetwork({
    
    Sys.sleep(2);
  
    datasetInput()
    
    
    
    
    
    
    
    
    
    
    
  })
  
  output$simple<-renderSimpleNetwork({
    
    layoutInput()
    
  })


  
  observe({
    visNetworkProxy("plot") %>%
      
   
      visOptions(highlightNearest = TRUE,nodesIdSelection=list(enabled=TRUE,values=input$nodeSelect))
  
   
    
    
  })
  observe({
    visNetworkProxy("plot") %>%
      
      visPhysics(enabled=input$phyactive)
    
    
  })
  observe({
    visNetworkProxy("plot") %>%
      
      visNodes(size=input$nodeSize)
    
    
  })
  
  
 observe({
    visNetworkProxy("plot") %>%
     
     visPhysics(solver=input$solve,barnesHut=list(gravitationalConstant=input$grav,centralGravity=input$centralGrav),repulsion=list(centralGravity=input$centralGrav),forceAtlas2Based=list(gravitationalConstant=input$grav,centralGravity=input$centralGrav))
   
   
      
  })


 

  ##VISNETWORK INPUT END##
  
  plotInput <- function(network,nodes,edges){
    
    cyNetwork <- network
    
    rashidcytoscapejs(nodes,edges,showPanzoom=TRUE)
  }
  
  
  
  output$downloadNetworkImage <- downloadHandler(
    filename = "Shinynet.csv",
    
    content = function(file) {
      statNet<-statNet[-c(1,3,4)]
      colnames(statNet)<-c("geneName","NCBI Information","geneID","Diff_Exp","score")
      statNet<-statNet[c("geneName","geneID","Diff_Exp","score","NCBI Information")]
      write.csv(statNet, file)
      
      
    })  
  output$downloadEnrich <- downloadHandler(
    filename = "Shinynet.csv",
    
    content = function(file) {
      
      write.csv(tableEnrichment, file)
      
      
    })  
  
  
  
  observeEvent(input$compDataSave, {
    filename=paste("files/",input$compData)
    write.csv(statNet,file=paste(filename,".csv",sep=""),row.names=FALSE,sep = "\t")
   
    session$sendCustomMessage(type = 'testmessage',
                              message = 'Thank you for clicking')
  })
  
  sig_name<-reactive({
   
    return(listfiles())
  })
  observeEvent(input$choices,{
    sigFileName<-paste("files/",input$choices,sep = "")
    sigFile<<-read.csv(file=sigFileName,sep='\t')
   
    
    
    
  })
  output$Signatures <- renderUI({
    selectInput("choices","Select Signatures",choices = as.vector(sig_name()))
  })
  output$text <- renderPrint({
    paste("you selected",input$choices,sep=" ")
  })
  
  
  
  
  
}
)
##other functions##
#list all files
listfiles <- function(){
  
  files <- Sys.glob("files/*.csv")
  
  l = list()
  for(i in 1:length(files)){
    
    l[[i]]=strsplit(files[i],"/")

  }
  
  vect = c()
  for(i in 1:length(files)){
   
    vect[i] = l[[i]][[1]][2]     
  }
  
  
  return (vect)
}
