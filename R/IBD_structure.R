#' identify the IBD structure of an SG (or Sliding-SG)
#'
#' @param gr.SG Grange of concerned SG
#' @param IBDClusters Dash output in data frame format
#' @param overlap the cutoff point at which the cluster is regarded as part of the structure
#'
#' @return Clusters of IBD structure
#' @export
#'
#' @examples IBD_structure(gr.SG, IBDClusters)
#' @examples IBD_structure(gr.SG=gr.SG, IBDClusters=df, overlap=0.5)
#' @examples IBD_structure(gr.SG=gr.Sliding_SG, IBDClusters=df, overlap=0.5)
#'
IBD_structure = function(gr.SG, IBDClusters, overlap = 0.5){
  chr = levels(seqnames(gr.SG))
  IBDClusters = IBDClusters[order(IBDClusters[,1]),   ]

  if (nrow(IBDClusters) > 0) {
    gr.IBDClusters = GRanges(seqnames = chr, ranges = IRanges(start = as.numeric(IBDClusters[,2]),
                                                              end = as.numeric(IBDClusters[,3])), score = score(gr.SG))
    ov = findOverlaps(gr.SG, gr.IBDClusters)
    gr.candidat = gr.IBDClusters[subjectHits(ov)]

    overlaps <- pintersect(gr.SG[queryHits(ov)], gr.candidat)
    percentOverlap <- round((width(overlaps)/width(gr.SG)), 2)

    IBD_structure.candidat = IBDClusters[,1][subjectHits(ov)]
    IBD_structure = IBD_structure.candidat[percentOverlap >  overlap]
    return(IBD_structure)
  }
}
