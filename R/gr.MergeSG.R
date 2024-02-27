#' Merging two SGs Grange
#'
#' @param gr.SG1 Grange of the first SG
#' @param gr.SG2 Grange of the second SG
#'
#' @return Merged Grange
#' @export
#'
#' @examples MergeSG(gr.SG1,gr.SG2)
#' @examples MergeSG(gr.Sliding_SG,gr.SG)
#'
#'
gr.MergeSG=function(gr.SG1, gr.SG2){
  start = min(start(gr.SG1), start(gr.SG2))
  end = max(end(gr.SG2), end(gr.SG1))
  if (levels(seqnames(gr.SG2)) == levels(seqnames(gr.SG1)) & score(gr.SG1)==score(gr.SG2)) {

    gr.SG1 = GRanges(seqnames = levels(seqnames(gr.SG1)),
                     ranges = IRanges(start = start, end = end),
                     score = score(gr.SG1))

    return(gr.SG1)
  }
}

