


shinyServer(function(input,output,session){
  
  
  ##HIDE/SHOW BASED ON WHETHER FILE UPLOADED OR EXAMPLE FILE LOADED
  observe({
       
      if ( input$load_example[1]==0 && is.null(input$file1) && session$clientData$url_search=="") {
        
        
        
        shinyjs::disable("downloadNetworkData")
              shinyjs::hide("saveImage")
              shinyjs::disable("downloadEnrich")
              shinyjs::hide("algorithm")
              shinyjs::hide("layout")
              shinyjs::hide("PPI")
              shinyjs::hide("goAn")
                 toggle(condition = FALSE  , selector = "#tabs li a[data-value=cytonet]")
                 toggle(condition = FALSE , selector = "#tabs li a[data-value=datTable]")
                 toggle(condition = FALSE , selector = "#tabs li a[data-value=irich]")
        
        
      }
    
      else{
        shinyjs::enable("downloadNetworkData")
        shinyjs::show("saveImage")
        shinyjs::enable("downloadEnrich")
        shinyjs::show("algorithm")
        shinyjs::show("layout")
        shinyjs::show("PPI")
        shinyjs::show("goAn")
        toggle(condition = TRUE  , selector = "#tabs li a[data-value=cytonet]")
        toggle(condition = TRUE , selector = "#tabs li a[data-value=datTable]")
        toggle(condition = TRUE , selector = "#tabs li a[data-value=irich]")
        # Hide the loading message when the rest of the server function has executed
               hide(id = "loading-content", anim = TRUE, animType = "fade")
        
      }
    
    })
  
  ##END ......HIDE/SHOW BASED ON WHETHER FILE UPLOADED OR EXAMPLE FILE LOADED
  
  
  
  
  
  
  ##LOAD FILE BASED ON WHETHER IT'S UPLOADED,LOAD EXAMPLE OR URL##
  
  observeEvent(session$clientData$url_search,{
    if(session$clientData$url_search!=""){
      
      nodeGO<<-NULL #To go back to default network construction and delete GO nodes
      edgeGO<<-NULL
      
      updateTabsetPanel(session, "tabs", selected = "cytonet")
    parsedGET<-parseQueryString(session$clientData$url_search)
     if (!is.null(parsedGET[["File"]])) {
        #path<-paste("http://www.ilincs.org/tmp/",parsedGET[["File"]],sep="")
       path<-paste("/mnt/raid/tmp/",parsedGET[["File"]],sep="")
        #path<-paste("/Users/Rashid/Desktop/Rashid/Career/PhD/Research/Events/BD2KAllHandsMeeting/signatures/",parsedGET[["File"]],sep="")
        loadedFile<<-read.csv(file=path,sep="\t")
      
        loadedFile<<-loadedFile[,1:5]
        colnames(loadedFile)<<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
       
    
    }
    
    }
    
  })
  
  observeEvent(input$load_example,{
    nodeGO<<-NULL #To go back to default network construction and delete GO nodes
    edgeGO<<-NULL
    updateTabsetPanel(session, "tabs", selected = "cytonet")
    loadedFile<<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
    loadedFile<<-loadedFile[,1:5]
    colnames(loadedFile)<<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
  })
  
  
  observeEvent(input$file1,{
    nodeGO<<-NULL #To go back to default network construction and delete GO nodes
    edgeGO<<-NULL
    updateTabsetPanel(session, "tabs", selected = "cytonet")
    loadedFile<<-read.csv(file=input$file1$datapath,sep="\t")
    loadedFile<<-loadedFile[,1:5]
    colnames(loadedFile)<<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
    
  })
  
  
  ##END.........LOAD FILE BASED ON WHETHER IT'S UPLOADED,LOAD EXAMPLE OR URL##
  
  
  ##GENERATE NETWORK PLOT##
  
  output$general_ui <- renderUI({
    
    visNetworkOutput("plot")
   
    
    
  })
  output$plot<-renderVisNetwork({
    
    datasetInput()
  })
  datasetInput<-reactive({
    
    algorithm<-input$algorithm
    interactionNetwork<-input$PPI
    networkLayout<-input$layout
    physics<-input$phyactive
    input$file1
    input$load_example
    input$goAn
    if(algorithm==1){
      
      ret<-RWR(loadedFile,phy=physics,layOut=networkLayout,package=FALSE,nodeGoData=nodeGO,edgeGoData=edgeGO,proteinN=interactionNetwork)
      
      netcon<<-ret$edgeData
      
      nodecon<<-ret$nodeData
      netData<<-ret$networkData #capture network data from algorithm
      
      ret$networkObj #return visnetwork object
     
    }
    else if(algorithm==3){
      
      ret<-dmGWAS(loadedFile,phy=physics,layOut=networkLayout,package=FALSE,nodeGoData=nodeGO,edgeGoData=edgeGO,proteinN=interactionNetwork)
      
      netcon<<-ret$edgeData
      
      nodecon<<-ret$nodeData
      netData<<-ret$networkData #capture network data from algorithm
      
      ret$networkObj #return visnetwork object
      
    }
    else if(algorithm==4){
      
      ret<-modifiedBioNetwork(loadedFile,phy=physics,layOut=networkLayout,package=FALSE,nodeGoData=nodeGO,edgeGoData=edgeGO,proteinN=interactionNetwork)
      
      netcon<<-ret$edgeData
      
      nodecon<<-ret$nodeData
      netData<<-ret$networkData #capture network data from algorithm
      
      ret$networkObj #return visnetwork object
      
    }
    
    else if(algorithm==5){
      
      ret<-topHundredNetwork(loadedFile,phy=physics,layOut=networkLayout,package=FALSE,nodeGoData=nodeGO,edgeGoData=edgeGO,proteinN=interactionNetwork)
      
      netcon<<-ret$edgeData
      
      nodecon<<-ret$nodeData
      netData<<-ret$networkData #capture network data from algorithm
      
      ret$networkObj #return visnetwork object
      
    }
    
  })
  
  
  ##END......GENERATE NETWORK PLOT##
  
  
  ##GET NETWORK DATA#
  output$contents<-renderDataTable({
    
    input$file1   # to have active reactive context everytime file changes
    input$load_example
    input$PPI
    input$algorithm
    
    netData<-netData[-c(1,3,4)]
    
    colnames(netData)<-c("geneName","NCBI_Information","geneID","Diff_Exp")
    netData<-netData[c("geneName","geneID","Diff_Exp","NCBI_Information")]
    
    
    netData$NCBI_Information<-paste0("<a href='",netData$NCBI_Information,"' target='_blank'>",netData$NCBI_Information,"</a>")
    netData # returning global variable storing network data 
    
    
    
  },escape=FALSE)
  
  
  ##END....GET NETWORK DATA##
  
  
  ##GET ENRICHMENT DATA##
  
  output$enrichmentAnalysis<-renderDataTable({
          input$file1   # to have active reactive context everytime file changes
          input$load_example
          input$PPI
          input$algorithm
          
          
          fdr<-input$enrichFDR
          
          
          #geninfo<-geneInfoFromPortals(geneList=as.character(loadedFile$GeneID),symbol=T,names=F)
          geninfo<- geneData[which(geneData$GeneID%in%as.character(loadedFile$GeneID)),]
          geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep="")) 
         print("enrichment")
      
           print(geneLabels)
          
          subnetGenes<-netData$geneID
          allGenes<-gsub("[\\(\\)]", "", regmatches(geneLabels, gregexpr("\\(.*?\\)", geneLabels)))
          
          result<-geneListEnrichment(subnetGenes, allGenes, functionalCategories = "GO", species = "Hs", minGenesInCategory = 10, maxGenesInCategory =1000, inBkg=TRUE, sigFDR = fdr , verbose=TRUE)
          
          
         
          
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
    
    
    
  })
  
  
  ##END...GET ENRICHMENT DATA##
  
  
  ##DOWNLOADS##
output$downloadNetworkData <- downloadHandler(
      filename = "",
      
      content = function(file) {
    
      write.csv(netData, file)
        
        
      })  
output$downloadEnrich <- downloadHandler(
        filename = "",
        
        content = function(file) {
        
        write.csv(tableEnrichment, file)

        
      })  
  ##END...DOWNLOADS##
  
 
  
  
  ##PHYSICS LAYOUT OPTIONS##
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
     
     
     ##END....PHYSICS LAYOUT OPTIONS##
     
     
##ADDING GO PATHWAYS TO NETWORK##
observeEvent(input$goAn,{

if(length(input$enrichmentAnalysis_rows_selected)>0){
idSel<-as.character(tableEnrichment[input$enrichmentAnalysis_rows_selected,1])
nameSel<-as.character(tableEnrichment[input$enrichmentAnalysis_rows_selected,2])
enrichment<-list("enrichID"=idSel,"enrichName"=nameSel)


enrichIDlength<-length(enrichment$enrichID)


from<-c()
to<-c()
for(i in 1:enrichIDlength){
addGO<-genesInEnrichedCategories(enrichment$enrichID[i], netData$geneID, funcCategories = "GO", species = "Hs")



for(j in 1:ncol(addGO)){
  
  if(addGO[,j]==TRUE){
    geneConversion<-substring( colnames(addGO)[j],2)
    GOnodeName<-geneConversion
    from<-c(from,enrichment$enrichID[i])
    to<-c(to,GOnodeName)
    
  
    
  
  
   
  }
  
}



}
  
  
if(input$PPI==1){
  GOdataframe<-data.frame(from,to,title=NA,weight=NA,X=NA,f.neighborhood=NA,f.fusion=NA,f.cooccurence=NA,f.coexpression=NA,f.experimental=NA,f.database=NA,f.textmining=NA,f.combined_score=NA,a_symbols=NA,b_symbols=NA) 
 
  }
else if (input$PPI==2){
GOdataframe<-data.frame(from,to) 
}


netconGO<-rbind(netcon,GOdataframe) #add from,to edgelist from existing network to newly added GO pathway

l<-nrow(nodecon)

for(i in 1:enrichIDlength){ ##add GO node attributes
  nodecon[l+i,1]<-as.character(enrichment$enrichID[i])
  nodecon[l+i,2]<-as.character(enrichment$enrichID[i])
  nodecon[l+i,3]<-"#ffff00 "
  nodecon[l+i,4]<-2
  nodecon[l+i,5]<-"black"
  nodecon[l+i,6]<-paste0("<p>",enrichment$enrichName[i],"</p>")
  nodecon[l+i,7]<-"square"
  nodecon[l+i,8]<-20
  
}





nodeGO<<-nodecon
edgeGO<<-netconGO

  



}
  
  else{
     #Remove GO nodes when no rows are selected
    nodeGO<<-NULL
    edgeGO<<-NULL
  }

})

##END...ADDING GO PATHWAY TO NETWORK

##DETECT CHANGES IN APP FOR REMOVING GO PATHWAYS##
    
observeEvent(input$algorithm,{
  
  nodeGO<<-NULL
  edgeGO<<-NULL
})

observeEvent(input$PPI,{
  
  nodeGO<<-NULL
  edgeGO<<-NULL
})

    
     
##END..DETECT CHANGES IN APP FOR REMOVING GO PATHWAYS##


##CHANGE TAB WHEN GO PATHWAYS ADDED/DELETED##

observeEvent(input$goAn,{
  

  updateTabsetPanel(session, "tabs", selected = "cytonet")
  
})
  
  
  
  

  
  
  
})


##END SERVER FUNCTION##





















# 
# library(DT)
# library(shiny)
# library(shinyjs)
# library(DLBCL)
# library(BioNet)
# library(igraph)
# library(CLEAN)
# library(CLEAN.Hs)
# library(networkD3)
# library(SigNetA)
# library(visNetwork)
# library(dmGWAS)
# #library(RCurl)
# 
# ##hello world!!!
# data(interactome)
# #statNet <- reactiveValues(df_data = NULL)
# #statNet<-statNet$df_data
# statNet<<-NULL
# visImg<<-NULL
# exFile<<-NULL
# en<<-NULL
# tableEnrichment<<-NULL
# shinyServer(function(input,output,session){
#   observe({
#     query<-parseQueryString(session$clientData$url_search)
#     
#     if (!is.null(query[["File"]])) {
#       File<<-query[['File']]
#       ##print(File)
#       #print("woohoo")
#       # File<<- read.delim(file=File)
#       # File<<- read.csv(file=File,sep="\t")
#       # #print(class(File))
#       # #print(data[[GeneNames]])
#       path<-paste("http://www.ilincs.org/tmp/",query[["File"]],sep="")
#      # path<-paste("http://dev.ilincs.org/tmp/",query[["File"]],sep="") #get file from url look it up
#      # path<-paste("/Users/Rashid/Desktop/Rashid/Career/PhD/Research/Events/BD2KAllHandsMeeting/signatures/",query[["File"]],sep="")
#       File<<-path
#       # File<<- read.csv(file=path,sep="\t")
#   #File<<-File[complete.cases(File),]
#   
#   shinyjs::enable("downloadNetworkImage")
#   shinyjs::show("saveImage")
#   shinyjs::enable("downloadEnrich")
#   shinyjs::show("algorithm")
#   shinyjs::show("layout")
#   shinyjs::show("PPI")
#   shinyjs::show("goAn")
#   toggle(condition = File  , selector = "#tabs li a[data-value=cytonet]")
#   toggle(condition = File, selector = "#tabs li a[data-value=datTable]")
#   toggle(condition = File, selector = "#tabs li a[data-value=irich]")
#   
# Sys.sleep(2)
#   
#   # Hide the loading message when the rest of the server function has executed
#   hide(id = "loading-content", anim = TRUE, animType = "fade")
#       ##print(path)
#       #File<<- read.delim(file=path)
#     }
#     
#     else{
#       File<<-""
#     }
#   })
#  observeEvent(input$load_example,{
#     
#     exFile<<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
#   
#    shinyjs::enable("downloadNetworkImage")
#    shinyjs::show("saveImage")
#    shinyjs::enable("downloadEnrich")
#    shinyjs::show("algorithm")
#    shinyjs::show("layout")
#    shinyjs::show("PPI")
#    shinyjs::show("goAn")
#    toggle(condition = TRUE  , selector = "#tabs li a[data-value=cytonet]")
#    toggle(condition = TRUE , selector = "#tabs li a[data-value=datTable]")
#    toggle(condition = TRUE , selector = "#tabs li a[data-value=irich]")
#    
#   })
#  #shinyjs::toggle(selector = "span.logo",TRUE,condition="[data-collapsed=true]" )
# 
#   #shinyjs::onclick("[role=button]",shinyjs::hide(selector = "[role=button]"))
#  
#  
#   observe({
#   #shinyjs::onclick(".skin-yellow .main-header .navbar .sidebar-toggle ",shinyjs::hide(".skin-yellow .main-header .logo "))
#     if (is.null(input$file1) || input$file1 == ""  ) {
#       
#         if(!(grepl("xls",File))){
#           
#       shinyjs::disable("downloadNetworkImage")
#       shinyjs::hide("saveImage")
#       shinyjs::disable("downloadEnrich")
#       shinyjs::hide("algorithm")
#       shinyjs::hide("layout")
#       shinyjs::hide("PPI")
#       shinyjs::hide("goAn")
#         }
#       else{
#         
#       }
#       
#     } else {
# 
#       shinyjs::enable("downloadNetworkImage")
#       shinyjs::show("saveImage")
#       shinyjs::enable("downloadEnrich")
#       shinyjs::show("algorithm")
#       shinyjs::show("layout")
#       shinyjs::show("PPI")
#       shinyjs::show("goAn")
#       
#       
#       Sys.sleep(2)
#       
#       # Hide the loading message when the rest of the server function has executed
#       hide(id = "loading-content", anim = TRUE, animType = "fade")
#       
#     }
#     
#     if(is.null(File) || File==""){
#     toggle(condition = input$file1  , selector = "#tabs li a[data-value=cytonet]")
#     toggle(condition = input$file1, selector = "#tabs li a[data-value=datTable]")
#     toggle(condition = input$file1, selector = "#tabs li a[data-value=irich]")
#     }
#     
#     
#     
#     
#   })
#   observeEvent(input$saveImage, {
#     
#     session$sendCustomMessage(type="saveImage", message="NULL")
#     
#   })
#  
# observe({
#   if(input$algorithm==3 || input$algorithm==4){
#     showModal(modalDialog(
#       title = "Important message",
#       "Please wait.This algorithm takes time to generate the subnetwork."
#     ))
#   }
# })
#   
#   output$clickedNode = renderPrint({
#     input$clickedNode
#   })
#   
#   output$connectedNodes = renderPrint({
#     input$connectedNodes
#   })
#   
#   ##gdgsgds
#   output$enrichmentAnalysis<-DT::renderDataTable({
#     if (is.null(input$file1) || input$file1 == "" ) {
#       if(!is.null(exFile)){
#         
#         enrichmentChange()
#         
#       }
#       else if(grepl("xls",File)){
#     enrichmentChange()
#      }
#       else{
#         paste("No genes analyzed")
#       }
#         
#       
#     } else if(!is.null(input$file) || !is.null(input$PPI)){
#       
#       
#      # if(class(err) == "try-error"){ 
#         
#        # output$input_error2=renderText("No significant enrichment analysis data generated.Please upload another signature.")
#       #  stopifnot(class(err) == "try-error")
#       #}
#      # else{
#        # output$input_error2=renderText("")
#      
#         enrichmentChange()
#         
#        #this is where enrichment functionality was
#         
#         
#       }
#       
#       
#   #  }
#   }
#   
#   
#   )
#   
#   #enrichmentChange<-reactive({
#    enrichmentChange<-function(){
#    # if((input$algorithm=="1" && input$layout=="1") || (input$algorithm=="1" && input$layout=="2"))
#     if(input$algorithm=="1")
#     {
#       if(is.null(input$file1)){
#         if(!(grepl("xls",File))){
#         logic<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
#         }
#         else{
#           logic<-read.csv(file=File,sep='\t')
#           colnames(logic)<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
#         }
#       }
#      
#       else{
#       logic<-read.csv(input$file1$datapath,sep="\t")
#       }
#       
#       
#       logic<-sortNetwork(logic)
#    
#       
#       geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
#       geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
#      
#       if(input$PPI=="2"){
#         subnetGenes<-statNet$id
#       }
#       else{
#         ##IMPORTANT CHANGES TO BE MADE HERE
#         #subnetGenes<-statNet$id ; Look into applying this
#       subnetGenes<-gsub("[\\(\\)]", "", regmatches( statNet$id, gregexpr("\\(.*?\\)",  statNet$id)))
#       }
#       allGenes<-gsub("[\\(\\)]", "", regmatches(geneLabels, gregexpr("\\(.*?\\)", geneLabels)))
#    
# 
# 
# 
#       result<-geneListEnrichment(subnetGenes, allGenes, functionalCategories = "GO", species = "Hs", minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, verbose=TRUE)
#     
#       
#    
#       
#       options(scipen = 999)
#       ID<-result$categories
#       Name<-result$Description
#       PValue<-format(result$FisherPValue,digits=2,scientific = TRUE)
#       FDR<-format(result$FisherFDR,digits=2,scientific = TRUE)
#       NetworkGenes<-result$nGenesInCategory
#       AllGenes<-result$nAllGenesInCategory
#       logOR<-round(result$logOR,2)
#       tableEnrichment<<-data.frame(ID,Name,PValue,FDR,NetworkGenes,AllGenes,logOR)
#      datatable(tableEnrichment, options=list(columnDefs=list(list(targets=4)),pageLength = 7))
#     
#       
#    }
#     
#    # else if((input$algorithm=="2" && input$layout=="2")||(input$algorithm=="2" && input$layout=="1")){
#     else if(input$algorithm=="2" || input$algorithm=="4"){
#       if(is.null(input$file1)){
#         if(!(grepl("xls",File))){
#           logic<-read.csv(file=system.file("extdata", "sig_try3.tsv", package = "SigNetA"),sep='\t')
#         }
#         else{
#           logic<-read.csv(file=File,sep='\t')
#           colnames(logic)<-c("signatureID","GeneID","GeneNames","coefficients","Pvals")
#         }
#       }
#       else{
#         logic<-read.csv(input$file1$datapath,sep="\t")
#       }
#       
#       
#       
#    
#       
#       geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
#       geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
#       
#       subnetGenes<-gsub("[\\(\\)]", "", regmatches( statNet$id, gregexpr("\\(.*?\\)",  statNet$id)))
#       allGenes<-gsub("[\\(\\)]", "", regmatches(geneLabels, gregexpr("\\(.*?\\)", geneLabels)))
#      
#       result<-geneListEnrichment(subnetGenes, allGenes, functionalCategories = "GO", species = "Hs", minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, verbose=TRUE)
#       
#    
#       
#       
#       options(scipen = 999)
#       ID<-result$categories
#       Name<-result$Description
#       PValue<-format(result$FisherPValue,digits=2,scientific = TRUE)
#       FDR<-format(result$FisherFDR,digits=2,scientific = TRUE)
#       NetworkGenes<-result$nGenesInCategory
#       AllGenes<-result$nAllGenesInCategory
#       logOR<-round(result$logOR,2)
#       tableEnrichment<<-data.frame(ID,Name,PValue,FDR,NetworkGenes,AllGenes,logOR)
#      
#       
#     }
#     
#    } 
#   #})
#   
#  # tableChange<-reactive({
#   tableChange<-function(){ 
# 
#    # if(input$algorithm=="1" ){
#     if(!is.null(statNet)){
#       
#       statNet<-statNet[-c(1,3,4)]
#       
#       colnames(statNet)<-c("geneName","NCBI Information","geneID","Diff_Exp")
#       statNet<-statNet[c("geneName","geneID","Diff_Exp","NCBI Information")]
#       
#       
#     }
#      
#       statNet
#    # }
#   }
#  # })
#   
#   output$contents <-DT::renderDataTable({
#     
#     if (is.null(input$file1) || input$file1 == "" ) {
#       if(!is.null(exFile) || is.null(input$PPI)){
#         
#         tableChange()
#         
#       }
#       else if(grepl("xls",File)){
#        
#           tableChange()
#       }
#       
#       else{
#         paste("No genes analyzed")
#       }
#     }
#     
#     else if(!is.null(input$file) || !is.null(input$PPI)){
#       
#       
#       
#       tableChange()
#       
#     }
#     
#     
#   })
# 
# 
#  observeEvent(input$algorithm,{
#   
#    tableChange()
#  })
#    
#  
#   observeEvent(input$goAn,{
#  
#     
#     idSel<<-as.character(tableEnrichment[input$enrichmentAnalysis_rows_selected,1])
#     nameSel<<-as.character(tableEnrichment[input$enrichmentAnalysis_rows_selected,2])
#     en<<-data.frame(idSel,nameSel)
#     
#   
#   if(nrow(en)>0){
#     
#    }
#    else{
#       
#       en<<-NULL
#      
#     }
#     
#     
#   
# #     output$plot<-renderVisNetwork({
# #     
# #       ret<-topHundredNetwork(input$file1$datapath,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=en,package=FALSE)
# #       
# #     })
# #    
#   
#   
#   })
#   output$genenodes <- renderText({ 
#     
#     if (is.null(input$file1) || input$file1 == "") {
#       paste("No genes analyzed")
#     } else {
#       
#       
#       newlist<-c()
#       
#       
#       for (i in 1:length(statNet$id)){
#         append(statNet$id,"\n",i)
#         newlist[i]<-statNet$id[i]
#       }
#       
#       newlist
#       
#       
#       
#       
#     }
#     
#     
#     
#     
#   })
#   
#   
#  
#    
#     
#     
# 
#   
#   output$general_ui <- renderUI({
#    if(input$algorithm=="1" && input$layout=="1"){
#      #rcytoscapejsOutput("plot", height="800px",width="900px")
#      visNetworkOutput("plot")
#    }
#     else if(input$algorithm=="2" &&  input$layout=="1"){
#       #rcytoscapejsOutput("plot", height="800px",width="900px")
#       visNetworkOutput("plot")
#     }
#     
#     
#     #LAYOUT
#      else if(input$algorithm=="1" && input$layout=="2"){
#       simpleNetworkOutput("simple")
#      }
#     else if(input$algorithm=="2" && input$layout=="2"){
#       simpleNetworkOutput("simple")
#     }
#     
#     else{
#       #rcytoscapejsOutput("plot", height="800px",width="900px")
#       visNetworkOutput("plot")
#     }
#     
#     
#     
#     
#   })
#   
#   enrichAdd<-function(){
#     
#     if(is.null(en)){
#        retu<-NULL;
#     }
#     else{
#       retu<-en
#     }
#   }
#     
#   
#   datasetInput <- reactive({
#     if(input$algorithm=="1" && (!grepl("xls",File))){
#       #if(!is.null(File))
#       #{ print(File)
#      #   ret<-topHundredNetwork(File,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=enrichAdd())
#       #}
#      # else{
#        
#      # print("algo 1")
#       RWR(input$file1$datapath,upload5="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
#    
#       #}
#       
#     }
#     else if(input$algorithm=="1" && (grepl("xls",File))){
#       #print("algo 1 1")
#       RWR(File,upload5="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
# 
#     }
#    
#     
#     
#     else if(input$algorithm=="2" && (!grepl("xls",File))){
#      # print("algo 2")
#        print(File)
#       ret<-bioNetwork(input$file1$datapath,upload2="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
#     }
#     
#     else if(input$algorithm=="2" && (grepl("xls",File))){
#      # print("algo 2 2")
#       ret<-bioNetwork(File,upload2="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
#     }
#     
#     else if(input$algorithm=="3"  && (!grepl("xls",File))){
#      # print("algo 3")
#       ret<-dmGWAS(input$file1$datapath,upload3="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
#     }
#     else if(input$algorithm=="3"  && (grepl("xls",File))){
#       ret<-dmGWAS(File,upload3="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
#     }
#     
#     else if(input$algorithm=="4" && (!grepl("xls",File))){
#      # print("algo4")
#       modifiedBioNetwork(input$file1$datapath,upload4="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
#     }
#     else if(input$algorithm=="4" && (grepl("xls",File))){
#      # print("algo4 4")
#       modifiedBioNetwork(File,upload4="yes",phy=input$phyactive,layOut=input$layout,package=FALSE)
#     }
#     else if(input$algorithm=="5"  && (!grepl("xls",File))){
#       ret<-topHundredNetwork(input$file1$datapath,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=enrichAdd())
#    
#     }
#     else if(input$algorithm=="5"  && (grepl("xls",File))){
# 
#       ret<-topHundredNetwork(File,upload1="yes",layOut=input$layout,proteinN=input$PPI,phy=input$phyactive,enrich=enrichAdd())
#     }
#     
#    
#     
#     
#     
#     
#     
#   })
#   
#   layoutInput<-reactive({
#     
#     #if(input$layout=="2"){
#       if(input$algorithm=="1"){
#         #source("D3layout.R")
#         ret<-D3layout(input$file1$datapath,upload3="yes")
#       }
#       else if(input$algorithm=="2"){
#       #source("D3bioNetwork.R")
#       ret<-D3bioNetwork(input$file1$datapath,upload4="yes")
#       }
#     }
#     
#     
#  # }
#     
#   )
#   
# 
#     
#     output$plot<-renderVisNetwork({
#     
#     Sys.sleep(2);
#   
#     datasetInput()
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#   })
#   
#   output$simple<-renderSimpleNetwork({
#     
#     layoutInput()
#     
#   })
# 
# 
#   
#   observe({
#     visNetworkProxy("plot") %>%
#       
#    
#       visOptions(highlightNearest = TRUE,nodesIdSelection=list(enabled=TRUE,values=input$nodeSelect))
#   
#    
#     
#     
#   })
#   observe({
#     visNetworkProxy("plot") %>%
#       
#       visPhysics(enabled=input$phyactive)
#     
#     
#   })
#   observe({
#     visNetworkProxy("plot") %>%
#       
#       visNodes(size=input$nodeSize)
#     
#     
#   })
#   
#   
#  observe({
#     visNetworkProxy("plot") %>%
#      
#      visPhysics(solver=input$solve,barnesHut=list(gravitationalConstant=input$grav,centralGravity=input$centralGrav),repulsion=list(centralGravity=input$centralGrav),forceAtlas2Based=list(gravitationalConstant=input$grav,centralGravity=input$centralGrav))
#    
#    
#       
#   })
# 
# 
#  
# 
#   ##VISNETWORK INPUT END##
#   
#   plotInput <- function(network,nodes,edges){
#     
#     cyNetwork <- network
#     
#     rashidcytoscapejs(nodes,edges,showPanzoom=TRUE)
#   }
#   
#   
#   
#   output$downloadNetworkImage <- downloadHandler(
#     filename = "Shinynet.csv",
#     
#     content = function(file) {
#       statNet<-statNet[-c(1,3,4)]
#       colnames(statNet)<-c("geneName","NCBI Information","geneID","Diff_Exp","score")
#       statNet<-statNet[c("geneName","geneID","Diff_Exp","score","NCBI Information")]
#       write.csv(statNet, file)
#       
#       
#     })  
#   output$downloadEnrich <- downloadHandler(
#     filename = "Shinynet.csv",
#     
#     content = function(file) {
#       
#       write.csv(tableEnrichment, file)
#       
#       
#     })  
#   
#   
#   
#   observeEvent(input$compDataSave, {
#     filename=paste("files/",input$compData)
#     write.csv(statNet,file=paste(filename,".csv",sep=""),row.names=FALSE,sep = "\t")
#    
#     session$sendCustomMessage(type = 'testmessage',
#                               message = 'Thank you for clicking')
#   })
#   
#   sig_name<-reactive({
#    
#     return(listfiles())
#   })
#   observeEvent(input$choices,{
#     sigFileName<-paste("files/",input$choices,sep = "")
#     sigFile<<-read.csv(file=sigFileName,sep='\t')
#    
#     
#     
#     
#   })
#   output$Signatures <- renderUI({
#     selectInput("choices","Select Signatures",choices = as.vector(sig_name()))
#   })
#   output$text <- renderPrint({
#     paste("you selected",input$choices,sep=" ")
#   })
#   
#   
#   
#   
#   
# }
# )
# ##other functions##
# #list all files
# listfiles <- function(){
#   
#   files <- Sys.glob("files/*.csv")
#   
#   l = list()
#   for(i in 1:length(files)){
#     
#     l[[i]]=strsplit(files[i],"/")
# 
#   }
#   
#   vect = c()
#   for(i in 1:length(files)){
#    
#     vect[i] = l[[i]][[1]][2]     
#   }
#   
#   
#   return (vect)
# }
