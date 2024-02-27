#' Computing Sall for a unique SG with disjoint clusters
#'
#' @param SG : One SG example  SG=SGList[[1]]
#'
#' @return  The value of Sall for this SG, a numeric value
#'
#' @export
#'
#' @examples  Sall_SG(SG=SGList[[1]])
#'
#'

Sall_SG= function(SG) {#SG=SGList[[1]]
  cls=SG[[2]] #[[1]][[2]]
  if (length(cls)>=1) {
    # Get list of subject IDs
    subjects <- unique(substring(unlist(cls), 1, regexpr("\\.", unlist(cls)) - 1))
    # Create empty list to store results
    subject_clusters <- list()
    # Loop over subject IDs and find corresponding clusters = create vector of alleles for each subject
    subject_clusters=lapply(subjects, function(subject) {
      cluster_names <- names(cls)[sapply(cls, function(x) subject %in% substring(x, 1, regexpr("\\.", x) - 1))]
      subject_clusters[[subject]] <- cluster_names
    })
    names(subject_clusters)=subjects



    # compute all possible sets consisting of one allele from each affected individual
    allele_sets <- expand.grid(subject_clusters)

    # get number of unique founder alleles
    founder_alleles <- unique(unlist(allele_sets))
    # compute product of number of times each founder allele appears in each row of allele_sets
    # allele_freqs <- apply(allele_sets, 1, function(x) {
    #   freqs <- sapply(founder_alleles, function(y) sum(x == y))
    #   prod(factorial(freqs))
    # })

    allele_freqs <- apply(allele_sets, 1, function(x) {
      freqs <- tabulate(match(x, founder_alleles), nbins = length(founder_alleles))
      prod(factorial(freqs))
    })

    Sall <- sum(allele_freqs)
  } else {
    Sall=1
    }
  return(Sall)
}
