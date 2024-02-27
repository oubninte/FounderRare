#'  @title Creating a Grange of Sliding_SG
#'
#'
#'  @description  Creating a Grange of Sliding_SG associate to a specific SG (to be specified). Length of Sliding_SG is twice the overlapSS (in Kb). For creating a more specific Sliding_SG gr.Next_SG function  give more flexibility
#'
#' @param gr.SG :  Grange  of concerned SG
#' @param overlapSS : Overlap between SG and Sliding_SG in Kbp : default value is 10 Kbp
#'
#' @return  : Grange of Sliding_SG
#' @export
#'
#' @examples Sliding_SG(gr.SG)
#' @examples Sliding_SG(gr.SG=gr.SG, overlapSS=10)
#' @examples Sliding_SG(gr.SG, 10)
#'
gr.Sliding_SG=function(gr.SG, overlapSS=10){
  overlapSS=overlapSS*1000
  start=end(gr.SG)-overlapSS
  SSG=GRanges(

    seqnames = levels(seqnames(gr.SG)) ,
    ranges = IRanges(start = start,
                     width = 2*overlapSS),
    score = score(gr.SG))

  return(SSG)
}


