#'
#'
#' convert the non-disjoint clusters in the SGlist output of TBM into the same output with segregated clusters.
#'
#' @param SGList The TBM function's output or any list with the same structure
#' @param IBDClusters Data frame format of dash output
#' @param remove Argument to exclude from the IBD structure any SG that does not have any clusters
#' @param disjM To construct disjoint clusters, you have two methods available. If you want to merge overlapping clusters, set disjM to "Merge". Alternatively, if you prefer to select the cluster with the most shared haplotype when overlaps occur, set disjM to "Bigger". Choose the method that best suits your data analysis needs.
#'
#' @return List of SG with disjoints clusters in IBD structure
#' @export
#'
#' @examples TBMD(SGList, IBDClusters)
#' @examples TBMD(SGList, IBDClusters, disjM="Merge", remove=TRUE)
#' @examples TBMD(SGList, IBDClusters, disjM="Bigger", remove=TRUE)
#'
#'
#'
TBMD=function(SGList, IBDClusters, disjM="Merge", remove=TRUE){
  toRemove=c()
  SGListD= list()


  if (disjM=="Bigger") {


            if (length(SGList) != 0) {
              for (i in 1:length(SGList)) {
                if (length(SGList[[i]][[2]])==0) {
                  toRemove=c(toRemove, i)
                }else {
                  SGListD[[i]]=Disj_Clst_SG_Bigger(SGList[[i]], IBDClusters)
                }
              }
            }

            if (remove & length(toRemove)>0) {
              SGListD=SGListD[-toRemove]
            }


      } else {
                if (disjM != "Merge") {
                  stop("\n The disjM parameter is not correct \n")
                }

                    if (length(SGList) != 0) {
                      for (i in 1:length(SGList)) {
                        if (length(SGList[[i]][[2]])==0) {
                          toRemove=c(toRemove, i)
                        }else {
                          SGListD[[i]]=Disj_Clst_SG(SGList[[i]], IBDClusters)
                        }
                      }
                    }

                    if (remove & length(toRemove)>0) {
                      SGListD=SGListD[-toRemove]
                    }


  }


  return(SGListD)
}

