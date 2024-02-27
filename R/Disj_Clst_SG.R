#' @title  Convert an SG's IBD structure into a disjoint clusters
#'
#'@description Use graph theory to convert an SG's IBD structure into a disjoint clusters
#'
#' @param SG Concerned SG
#' @param IBDClusters The output of Dash, which is used to build SG, shows which clusters share the same haplotype (IBD)
#'
#' @return SG with Disjoint clusters
#' @export
#'
#' @examples Disj_Clst_SG(SG, IBDClusters)
#'
#'
#'
Disj_Clst_SG=function(SG, IBDClusters){
  columns_to_remove <- seq(4, ncol(IBDClusters), by = 2)
  IBDClusters <- IBDClusters[, -columns_to_remove]
  IBD_Struc=SG[[2]]

  # Subset the data frame and select relevant columns
  tmp <- IBDClusters[IBDClusters[, 1] %in% IBD_Struc, c(1, 4:ncol(IBDClusters))]

  # Remove columns with all NA values
  tmp2 <- tmp[, colSums(is.na(tmp)) != nrow(tmp)]

  # Transpose the data frame and set column names
  tmp2 <- as.data.frame(t(tmp2)[-1,])

  # Remove leading and trailing spaces from column names
  colnames(tmp2) <- gsub(" ", "", colnames(tmp2))

  # Remove NAs and leading/trailing spaces from data
  tmp3 <- lapply(tmp2, function(x) gsub(" ", "", as.vector(na.omit(x))))

  # Generate edges
  edges <- do.call(rbind, lapply(tmp3, function(x) {
    if (length(x) > 1) cbind(head(x, -1), tail(x, -1)) else NULL
  }))


  if (is.null(edges)) {
    # If edges are null, create a list with tmp2
    IBD_Struc_Disj <- list(cd1 = tmp2)
  } else {
    # If edges are not null, create a graph and split into clusters
    g <- graph.data.frame(edges, directed = FALSE)
    IBD_Struc_Disj <- split(V(g)$name, clusters(g)$membership)
    names(IBD_Struc_Disj) <- paste0("cd", seq_along(IBD_Struc_Disj))
  }


  SG[[2]]=IBD_Struc_Disj

  return(SG)
}
