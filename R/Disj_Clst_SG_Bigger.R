#' @title  Convert an SG's IBD structure into a disjoint clusters
#'
#' @description This code first sorts the groups by size in descending order. Then it checks each group for overlap with the groups already. When there are overlaps, the largest group is kept.
#'
#' @param SG Concerned SG
#' @param IBDClusters The output of Dash, which is used to build SG, shows which clusters share the same haplotype (IBD)
#'
#' @return SG with Disjoint clusters
#' @export
#'
#' @examples Disj_Clst_SG_Bigger(SG, IBDClusters)
#'
#'
#'
Disj_Clst_SG_Bigger=function(SG, IBDClusters){
  columns_to_remove <- seq(4, ncol(IBDClusters), by = 2)
  IBDClusters <- IBDClusters[, -columns_to_remove]
  IBD_Struc=SG[[2]]
  IBDClusters=IBDClusters[IBDClusters[,1] %in%IBD_Struc,]
  IBD_Struc_list <- lapply(1:nrow(IBDClusters), function(i){
    # Get the name for the list element from the first column
    element_name <- as.character(IBDClusters[i, 1])
    # Get the elements from column 4 to the end, excluding NA values
    elements <- IBDClusters[i, 4:ncol(IBDClusters)]
    elements <- elements[!is.na(elements)]
    gsub(" ", "", elements)
  } )
  names(IBD_Struc_list)=IBD_Struc
  IBD_Struc_list <- IBD_Struc_list[order(sapply(IBD_Struc_list, length), decreasing = TRUE)]

  # Create a function to check for overlaps
  check_overlap <- function(current_group, existing_groups) {
    overlap <- FALSE
    for (existing_group in existing_groups) {
      if (length(intersect(current_group, existing_group)) > 0) {
        overlap <- TRUE
        break
      }
    }
    return(!overlap)
  }

  # Initialize an empty list to hold the non-overlapping groups
  non_overlapping_groups <- list()

  # Iterate over each group
  for (current_group in IBD_Struc_list) {
    # If the current group does not overlap with any group in the non_overlapping_groups list, add it to the list
    if (check_overlap(current_group, non_overlapping_groups)) {
      non_overlapping_groups <- c(non_overlapping_groups, list(current_group))
    }
  }

  names(non_overlapping_groups) <- paste0("cd", seq_along(non_overlapping_groups))

  SG[[2]]=non_overlapping_groups
  return(SG)


}
