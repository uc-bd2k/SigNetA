#' Sort Network Function
#'
#' This function allows you to sort top 100 genes from a signature 
#' @param File ? Defaults to NULL.
#' @keywords sortNetwork
#' @export
#' @examples
#' sortNetwork()  will sort network for any signature file
sortNetwork<-function(x){
#File<- read.csv(file=x,sep="\t")

sorted<-x[order(x$Pvals),] 

topgenes<-head(sorted,100)

}
