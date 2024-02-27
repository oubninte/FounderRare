#'@title  Function based on the tunnel-boring machine (TBM) principle to create one SG;
#'
#'@description  Function based on the tunnel-boring machine (TBM) principle to create one SG;
#'primarily designed to expand the existing SG if its IBD structure like that of a sliding-SG;
#' else it produces a new nearby SG
#'
#' @param SG Concerned SG
#' @param IBDClusters IBD clusters (formatted as a data frame from Dash output)
#' @param overlapSS Overlap between Sliding-SG and SG grange; 10Kbp is the default value
#' @param overlapNS Overlap between the SG grange and the Next-SG grange; by default, this is set to 0
#' @param overlapI the threshold at which the cluster is considered to be included in the IBD structure
#' @param MinL   The SG Grange's minimum length in Kbp is 400 Kbp by default
#'
#' @return A list of two elements; the Grange of SG and IBD structure
#' @export
#'
#' @examples TBM_SG(SG, IBDClusters)
#' @examples TBM_SG(SG, IBDClusters, overlapSS=10, overlapNS=0, overlapI=0.5, MinL=400)
#'
#'
#'
#'
TBM_SG=function(SG, IBDClusters, overlapSS = 10, overlapNS = 0, overlapI = 0.5, MinL = 400){

  gr.SG = SG[[1]]
  SG_IBDstructure = SG[[2]]
  gr.Sliding_SG = gr.Sliding_SG(gr.SG, overlapSS)
  SSG_IBDstructure = IBD_structure(gr.Sliding_SG, IBDClusters,   overlapI)

  Sliding_SG = list(gr.Sliding_SG, SSG_IBDstructure)
  Action=""

  if (is_same_IBDstructure(Sliding_SG, SG)) {
    gr.SGN = gr.MergeSG(gr.SG, gr.Sliding_SG)
    SGN_IBDstructure=SG_IBDstructure
    #SGN = list(gr.SGN, SG_IBDstructure)
    Action="Merge"
  }  else {
    gr.SGN = gr.Next_SG(gr.SG, MinL, overlapNS)
    SGN_IBDstructure = IBD_structure(gr.SGN, IBDClusters, overlapI)
    #SGN = list(gr.SGN, SGN_IBDstructure)
    Action="New"
  }

  SGN = list(gr.SGN, SGN_IBDstructure)
  return(c(SGN , Action))
}
