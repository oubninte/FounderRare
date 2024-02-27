#' Computing Sall for all SGs of chromosome
#'
#' @param SGList list of SG with disjoints IBD clusters for a specific chromosome
#' @param Naff Number of affected individuals
#'
#' @return Dataframe with the value of Sall for each SG, Sall proposed by Whittemore and Halpernand and the Sall_WH per affected if Naff is gived in input
#' @export
#'
#' @examples Sall(SGList)
#' @examples Sall(SGList[[1]])
#' @examples Sall(SGList[1:5])
#'
#'
#'
Sall=function(SGList, Naff=NULL){


  Df_Sall = data.frame(chr = character(), SG = character(), Start = double(), End = double(), Sall = double())

  Df_Sall=rbind(Df_Sall, t(sapply(1:length(SGList), function(k){
  data.frame(chr = levels(seqnames(SGList[[k]][[1]])),
  SG=paste0("SG", k),
  Start=start(SGList[[k]][[1]]),
  End=end(SGList[[k]][[1]]),
  Sall=Sall_SG(SGList[[k]])
  )  } )))

    if (!is.null(Naff)) {
    Df_Sall$Sall_WH= as.numeric(Df_Sall$Sall)/(2^Naff)
    Df_Sall$SallWH_perAff= as.numeric(Df_Sall$Sall_WH)/Naff
  }
  outSall=as.data.frame(unnest(Df_Sall, cols = colnames(Df_Sall)))
  return(outSall)
}

