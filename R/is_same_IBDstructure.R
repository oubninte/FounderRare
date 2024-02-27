#' Compare two SGs' IBD structures
#'
#' @param SG1 The first SG: IBD_structure is the second element in SG, which is a list of two elements.
#' @param SG2 The second SG
#'
#' @return Return True if SG1 and SG2 have the same IBD_structure otherwise False
#' @export
#'
#' @examples is_same_IBDstructure(SG1, SG2)
is_same_IBDstructure= function(SG1, SG2 ){

  IBD_structure1=SG1[[2]]
  IBD_structure2=SG2[[2]]
  return(identical(IBD_structure1[order(IBD_structure1)],IBD_structure2[order(IBD_structure2)]))
}
