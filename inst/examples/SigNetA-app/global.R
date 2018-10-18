library(DT)
library(shiny)
library(shinyjs)
library(DLBCL)
library(BioNet)
library(igraph)
library(CLEAN)
library(CLEAN.Hs)

library(SigNetA)
library(visNetwork)
#library(dmGWAS) #version 2.4
library(RCurl)
library(devtools)

#library(roxygen2)
#data(interactome)

#document()

loadedFile<<-c()  #global File input handler
netData<<-c() #global network data input handler
netcon<<-c() #connected from,to in network data
nodeGO<<-NULL #visobject node  with GO terms in it
edgeGO<<-NULL #visobject edge with GO terms in it

##Load STRING Data##


#load(url("http://www.ilincs.org/tmp/signeta/weightedGraphStringPPI_10.rda")) #STRING PPI network
#load(url("http://www.ilincs.org/tmp/signeta/stringToEntrezDetailed.Rdata")) #edge data

#For local use
#load( url("http://www.ilincs.org/tmp/signeta/lincscp_1.rda"))
# load("/Users/Rashid/Desktop/Rpackages/SigNetA_Data/weightedGraphStringPPI_10.rda")#new STRING PPI network
# load("/Users/Rashid/Desktop/Rpackages/SigNetA_Data/interactome.rda")
# load("/Users/Rashid/Desktop/Rpackages/SigNetA_Data/geneData.rda")

#For server use#

load("/tmp/datasets/weightedGraphStringPPI_10.rda") #STRING PPI network
#load("/tmp/datasets/stringToEntrezDetailed.rda") #edge data
load("/tmp/datasets/interactome.rda")#interactome
load( "/tmp/datasets/lincscp_1.rda")
load("/tmp/datasets/geneData.rda")



##END..Load STRING Data##



