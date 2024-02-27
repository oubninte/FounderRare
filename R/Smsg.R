#' Function to compute Smsg
#'
#' @param SGList List of SG with disjoints IBD clusters for a given chromosome
#' @param Naff Number of affected individuals
#'
#' @return Dataframe with the value of Smsg for each SG and, in the event that Naff is provided as input, the Smsg per affected also is calculated
#' @export
#'
#'
#'
#' @examples Smsg(SGList)
#' @examples Smsg(SGList[[1]])
#' @examples Smsg(SGList[1:5])
#'
#'
#'
Smsg=function(SGList, Naff=NULL){


  Df_Smsg = data.frame(chr = character(), SG = character(), cd = character(),
                       Start = double(), End = double(), Smsg = double())
  Df_Smsg = rbind(Df_Smsg, t(sapply(1:length(SGList),
                                    function(k) {
                                      SG_k = SGList[[k]][[2]]
                                      S_i = sapply(1:length(SG_k), function(cd) length(SG_k[[cd]]))
                                      cls_keep = which(S_i == max(S_i))
                                      data.frame(chr = levels(seqnames(SGList[[k]][[1]])),
                                                 SG = paste0("SG", k), cd = names(SG_k)[cls_keep],
                                                 Start = start(SGList[[k]][[1]]), End = end(SGList[[k]][[1]]),
                                                 Smsg = max(S_i))
                                    })))
  Df_Smsg <- unnest(Df_Smsg, cols = colnames(Df_Smsg))

  if (!is.null(Naff)) {
    Df_Smsg$Smsg_per_aff = as.numeric(Df_Smsg$Smsg)/Naff
  }

  return(as.data.frame(Df_Smsg))
}

