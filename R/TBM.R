#' Function based on the tunnel-boring machine (TBM) principle to creat a list of SGs
#'
#'  Based on TBM_SG function to scan all chromosome
#'
#' @param SGList list of SGs
#' @param IBDClusters IBD clusters (formatted as a data frame from Dash output)
#' @param overlapSS Overlap between Sliding-SG and SG grange; 10Kbp is the default value
#' @param overlapNS Overlap between the SG grange and the Next-SG grange; by default, this is set to 0
#' @param overlapI the threshold at which the cluster is considered to be included in the IBD structure
#' @param MinL   The SG Grange's minimum length in Kbp is 400 Kbp by default
#' @param N Number of iteration, The NULL value is used to scan the entire chromosome.
#' @param chr Concerned chromosome
#'
#' @return list of SGs
#' @export
#'
#' @examples TBM(SGList, IBDClusters)
#' @examples TBM(SGList, IBDClusters, overlapSS=10, overlapNS=0, overlapI=0.5, MinL=400)
#'
#'
TBM=function(IBDClusters, SGList=list(), chr=NULL, N=NULL, overlapSS = 10, overlapNS = 0, overlapI = 0.5,MinL = 400){
  if (length(SGList)==0 ) {

    if (is.null(chr)) {
      chr="*"
    }

    gr.SG1=GRanges(
      seqnames = chr ,
      ranges = IRanges(start = min(IBDClusters[,2]),
                       width = MinL*1000),
      score = "*")
    SG1= list(gr.SG1, IBD_structure(gr.SG1, IBDClusters))
    SGList=list(SG1)
  }


  if (is.null(N)) {
    while(end(SGList[[length(SGList)]][[1]])<max(IBDClusters[,3]) ) {
      SG=TBM_SG(SGList[[length(SGList)]], IBDClusters, overlapSS, overlapNS, overlapI, MinL)
      if (SG[[3]]=="New") {
        SGList[[length(SGList)+1]]=SG[1:2]
      } else {SGList[[length(SGList)]]=SG[1:2]}

    }
  } else{
    if (N>0) {
      for (i in 1:N ){
        SG=TBM_SG(SGList[[length(SGList)]], IBDClusters, overlapSS, overlapNS, overlapI, MinL)
        if (SG[[3]]=="New") {
          SGList[[length(SGList)+1]]=SG[1:2]
        } else {SGList[[length(SGList)]]=SG[1:2]}
      }
    }
  }

  return(SGList)
}
