#' Construct the Next SG range for a particular SG.
#'
#' @param gr.SG Grange of concerned SG
#' @param minL Minimum length of SG Grange in Kbp, 400Kbp by default
#' @param Overlap  Overlap between current SG and next SG in bp; equal to 0 by default
#'
#' @return Grange of next SG
#' @export
#'
#' @examples gr.Next_SG(gr.SG)
#' @examples gr.Next_SG(gr.SG, minL=400)
#' @examples gr.Next_SG(gr.SG, minL=400, Overlap=0)
#'
#'
gr.Next_SG=function (gr.SG, minL = 400, Overlap = 0) {
  start = end(gr.SG) - Overlap +1
  NSG = GRanges(seqnames = levels(seqnames(gr.SG)),
                ranges = IRanges(start = start, width = minL * 1000),
                score = score(gr.SG))
  return(NSG)
}
