# FounderRare: A Package for Identifying Rare Variants in Complex Diseases

## Introduction

**FounderRare** is a package designed to provide a deeper understanding of the novel Smsg method and its potential applications in the field of complex disease research. This comprehensive guide elucidates the workings of the package, expounds on the functions employed, and provides illustrative examples of the output.

## About Smsg

Smsg is a groundbreaking statistical approach that leverages Founder Population Genealogy and Identical-By-Descent (IBD) Segments to identify rare variants in complex diseases. This innovative method employs IBD segments as markers of recent variants within family data derived from a founder population with accessible genealogy.

## How It Works

Smsg ingeniously partitions the genome into synthetic genes (SGs), identifies clusters of affected individuals who share a specific IBD segment over an SG, and conducts statistical tests to assess the enrichment of SG sharing among affected individuals.

## Scalability and Simulation

Capitalizing on genotype array data, Smsg is scalable to large sample sizes. Furthermore, it harnesses msprime to simulate the transmission of the entire genome within a genealogy, thereby determining the null distribution of the statistics.

## Conclusion

With **FounderRare**, researchers and scientists can now leverage the power of Smsg to identify rare variants in complex diseases, opening up new avenues for understanding and treating these conditions.
