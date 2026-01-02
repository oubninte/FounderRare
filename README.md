# FounderRare: A Package for Identifying Rare Variants in Complex Diseases

## Introduction

**FounderRare** is a package designed to provide a deeper understanding of the novel Smsg method and its potential applications in the field of complex disease research. Smsg is a groundbreaking statistical approach that leverages Founder Population Genealogy and Identical-By-Descent (IBD) Segments to identify rare variants in complex diseases. This innovative method employs IBD segments as markers of recent variants within family data derived from a founder population with accessible genealogy.

## How It Works
/vignettes/About_FounderRare.Rmd is comprehensive guide elucidates the workings of the package, expounds on the functions employed, and provides illustrative examples of the output.
More information are availables :
Articles url ()

###
If you use this package please cite this article or this one

## Licence
[MIT Licence](LICENCE)

# Introduction

In this tutorial, we elucidate the workings of `FounderRare` Package, expound on the functions employed, and provide illustrative examples of the output. This comprehensive guide aims to provide a deeper understanding of the novel Smsg method and its potential applications in the field of complex disease research.

Smsg is a groundbreaking approach, Statistical approach leveraging Founder Population Genealogy and IBD Segments to Identify Rare Variants in Complex Diseases. This innovative method employs identical-by-descent segments (IBDS) as markers of recent variants within family data derived from a founder population with accessible genealogy.

Smsg ingeniously partitions the genome into synthetic genes (SGs), identifies clusters of affected individuals who share a specific IBD segment over an SG, and conducts statistical tests to assess the enrichment of SG sharing among affected individuals.

Capitalizing on genotype array data, Smsg is scalable to large sample sizes. Furthermore, it harnesses msprime to simulate the transmission of the entire genome within a genealogy, thereby determining the null distribution of the statistics.



## R Markdown

This is an R Markdown document for `FounderRare` package. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the `Knit` button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. 

## Downald and install FounderRare Package

The packages `devtools`, `BiocManager` and `GenomicRanges` are essential prerequisites for the installation of the `FounderRare` package. Additionally, the utilization of `FounderRare` necessitates the presence of `dplyr`, `igraph`, `stringr`, and `tidyr`. It is highly advised to install and load all these packages prior to attempting to use the `FounderRare` package.

The corresponding code below, demonstrating a method to install the `FounderRare` package:

    # This line is used to install the 'devtools' package. 'devtools' provides tools to make developing R packages easier.
     install.packages("devtools")
    
    # This block of code checks if the 'BiocManager' package is installed. If not, it installs the package. 'BiocManager' is needed to install Bioconductor packages.
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    # This line is used to install the 'GenomicRanges' package using 'BiocManager'. 'GenomicRanges' provides efficient data structures and methods for manipulating genomic intervals and variables defined along a genome.
    BiocManager::install("GenomicRanges")
    
    # This line is used to install the 'FounderRare' package from GitHub using 'devtools'. You need to replace "username" with the actual GitHub username where the package is hosted and adjust the path accordingly.
    devtools::install_github("oubninte/FounderRare") 


Again I recommend to correctly loaded the dependencies libraries, here a code for that:

Ensure that the placeholder `oubninte/` in the devtools::install_github() function invocation is substituted with the accurate GitHub username where the `FounderRare` package is located. It is important to note that if `FounderRare` is not hosted on GitHub or if the package structure deviates, it may be necessary to alter the installation command to suit the specific circumstances.

I advocate for the correct loading of the dependency libraries. The following code can be used for this purpose:
```{r echo=TRUE, message=FALSE, warning=FALSE}
# Load the 'GenomicRanges' package, which provides efficient data structures and methods for manipulating genomic intervals and variables defined along a genome.
library(GenomicRanges)

# Load the 'FounderRare' package. It seems you've loaded this package twice in your code. You only need to load it once.
library(FounderRare)

# Load the 'igraph' package, which provides tools for network analysis.
library(igraph)

# This is a duplicate line and is not necessary, as you've already loaded the 'FounderRare' package above.
library(FounderRare)

# Load the 'dplyr' package, which provides a flexible toolkit for data manipulation in R.
library(dplyr)

# Load the 'tidyr' package, which provides functions for tidying data.
library(tidyr)
```

## Enumeration of Functions Contained in the Package

A robust way to list all functions in a package is to use the ls() function along with the pattern argument. This code uses the ls() function with the pattern argument set to `\^[[:alpha:]]`, which matches names starting with a letter (typical for function names). The "package:FounderRare" specification ensures that it looks specifically within the "FounderRare" package.

Note: Make sure you have loaded the package using library(FounderRare) before running this code.

```{r}
# List all functions in the package
all_functions <- ls("package:FounderRare", pattern = "^[[:alpha:]]")

# Print the names of all functions
print(all_functions)
```

## Loading and importing data 

The FounderRare package incorporates the IBDClusters dataset as a key component, serving as a representative example of the data format required for the package’s input arguments. The dataset is loaded into the environment using the data("IBDClusters") command. To gain a comprehensive understanding of the dataset’s structure, the str(IBDClusters) command is utilized. Notably, this dataset is generated by parsing Dash’s .hcl file.

Here’s an analytical interpretation of the output:

    The dataset, named "IBDClusters," is a data frame with 27 rows and 13 columns.
    The columns represent various information related to IBD (Identity by Descent) clusters.
    The first column is the cluster identifier.
    Columns 2 and 3 represent the start and end positions of the clusters.
    Every alternate pair of columns (such as 4 and 5) convey the Family and Sample ID data for segment sharers.

Executing head(IBDClusters) will display the initial few rows of the dataset, thereby facilitating a preliminary inspection of the actual data.

```{r }
# Incorporation of data
data("IBDClusters")

# Overview of the dataset structure
str(IBDClusters)

# Preview of the dataset
head(IBDClusters)
```

**Importing data from output of Dash**

The specificity of importing file.hcl in this code lies in the dynamic generation of column names based on the number of affected individuals . The read.table function is used to read the file and store it as a data frame. The col.names argument is used to specify the column names of this data frame. In this case, the column names are generated by concatenating **V** with a sequence of numbers. The length of this sequence is **2*2*Naff+3**, this opération is needed to ensure that the file is imported correctly. The fill=T argument ensures that if the rows have unequal length, blank fields are filled with NA. This is particularly useful when dealing with data files that may have missing values.


    # Number of affected individuals
    Naff=19 
    
    # Path to the .hcl file
    clusterIBD_path="path/to/Dash/File_chr22.hcl" 
    
    # Read the .hcl file into a table. The column names are generated dynamically based on the number of affected individuals (Naff).
    # The 'fill' argument is set to TRUE to fill in missing values with NA.
    IBDClusters=read.table(clusterIBD_path, col.names = paste0("V",seq_len(2*2*Naff+3)), fill=T)



# Constructing Initial Genomic Ranges (Granges) for the Primary Synthetic Gene (SG) and Sliding-SG

It is always advisable to establish the initial Grange prior to utilizing the package function. This approach allows for the inclusion of certain unavailable data, such as the chromosome number. The following information is required to create the first element:

    chr = 22: Denotes the chromosome number (in this instance, chromosome 22).
    minL = 400: Indicates the minimum length of SG in kilobases (Kbp).
    GRanges(...): Constructs a genomic range object.
        seqnames = chr: Designates the chromosome for the range.
        ranges = IRanges(...): Defines the start and width of the range.
        score = "*": No metadata.


```{r }
chr=22
minL=400 
# Creation of a Grange object for the first synthetic gene (SG)
gr.SG1=GRanges(
  seqnames = chr , # Chromosome for the range
  # Start and width of the range are defined
  ranges = IRanges(start = min(IBDClusters[,2]),
                   width = minL*1000),
  score = "*") # No metadata


```

This code block creates a genomic range object for the first synthetic gene (SG). The GRanges function from the GenomicRanges package in R is used to create this object. The seqnames argument specifies the chromosome for the range (in this case, chromosome 22), and the ranges argument, which uses the IRanges function, specifies the start and width of the range. The score argument is set to "\*", indicating that there is no metadata associated with this genomic range. The start of the range is set to the minimum value of the second column of the IBDClusters data frame, and the width of the range is set to minL*1000, where minL is the initial length of SG in kilobases (Kbp).

**Generation of a GRange for the Sliding Synthetic Gene (gr.Sliding_SG):**

The gr.Sliding_SG function is specifically designed to create a GRange object, representing a sliding synthetic gene (Sliding-SG), which corresponds to a given synthetic gene GRange (gr.SG). As illustrated in the provided code snippet, the function is invoked with the parameters gr.SG and overlapSS = 10 (default value). In essence, this function generates a sliding synthetic gene (SSG) with a length that is double the specified overlap (overlapSS), measured in kilobases (Kbp). The overlapSS parameter offers researchers the ability to customize the window size as per their research requirements.


```{r }

# The gr.Sliding_SG function is called with parameters gr.SG1 and overlapSS = 10
# This generates a sliding synthetic gene (SSG) with a length that is twice the specified overlap (overlapSS)
gr.Sliding_SG(gr.SG1, overlapSS = 10)

```

# Generation of the Subsequent Genomic Range via *gr.Next_SG*

The `gr.Next_SG` function is designed to construct a Grange of subsequent synthetic gene (SG) following a specified SG Grange (gr.SG). This function provides customization options, enabling users to stipulate the minimum length (minL) of the ensuing SG in kilobases. Moreover, an optional parameter, Overlap, allows for the creation of an overlap between the current and subsequent SGs. The code sample demonstrates the utilization of this function by invoking it with parameters gr.SG, minL = 400, and Overlap = 0. This adaptability allows researchers to tailor the genomic segment analysis to meet specific research needs. 

Note: 
  This function can be employed to generate a custom sliding SG.
  In this instance, invoking gr.Next_SG(gr.SG1) will yield the same result.  minL = 400 and Overlap = 0 are default values

```{r }
# The gr.Next_SG function is called with parameters gr.SG1, MinL = 400, and Overlap = 0
# This generates a subsequent synthetic gene (SG) Grange with a specified minimum length (MinL) and an optional overlap
gr.Next_SG(gr.SG1, minL = 400, Overlap = 0)

```


# Analysis of Identity by Descent (IBD) Structure within a Synthetic Gene 

The provided R code presents a test case that involves the detection and comparison of the IBD structures within two SGs (exactly between SG and sliding-SG). 

##Detection of the IBD Structure within an SG (IBD_structure):

The `IBD_structure` function accepts a Grange (gr.SG), a data frame representing IBD clusters (IBDClusters), and an optional overlap cutoff that represents the percentage of overlap between the cluster and SG at which a cluster is considered to be included in the IBD Structure (default value is 0.5). The function returns the clusters of the IBD structure based on the specified parameters, leveraging information from the IBD clusters and calculating the percentage overlap with the given genomic range.

```{r }
# The IBD_structure function is called with parameters gr.SG1 and IBDClusters
# This generates the IBD structure of the first synthetic gene (SG1)
IBD_structure(gr.SG1, IBDClusters) 


# The gr.Sliding_SG function is called with parameter gr.SG1
# This generates a sliding synthetic gene (SSG) Grange
SG1= list(gr.SG1, IBD_structure(gr.SG1, IBDClusters))
gr.SSG=gr.Sliding_SG(gr.SG1)
SSG=list(gr.SSG, IBD_structure(gr.SSG, IBDClusters))

```

## Comparison of IBD Structures :
        The `is_same_IBDstructure` function compares the IBD structures of two SGs (SG1 and Sliding-SG in this example) to ascertain if they are identical. The function returns TRUE if the IBD structures are identical and FALSE otherwise.

```{r }
# The is_same_IBDstructure function is called with parameters SG1 and SSG
# This compares the IBD structures of SG1 and SSG
is_same_IBDstructure(SG1, SSG)

```

# Running the *TBM* Function

The R code presented below serves as a test case for the TBM (from Tunnel-Boring Machine) function. This function is instrumental in identifying and comparing Identity by Descent (IBD) structures. The TBM function is predicated on the TBM_SG function, which is used to generate a list of SG.

## Testing the *TBM_SG* Function:

The `TBM_SG` function, built on the TBM principle, dynamically creates or expands Granges based on their IBD structures. It accepts several parameters, including the SG under consideration, IBD clusters (IBDClusters), and various overlap and length parameters.

If the IBD structure of the existing SG mirrors that of a Sliding-SG, as determined by the is_same_IBDstructure function, the function augments the existing SG. Specifically, it generates a sliding-SG and compares its IBD structure with that of the input SG. If they align, the function merges the sliding SG with the original SG, returning a list containing the extended SG Grange and its IBD structure, with an action indicator set to “Merge”. If the IBD structures diverge, indicating a unique pattern, the function generates a new proximate Grange using the gr.Next_SG function. This new SG is predicated on the minimum length parameter (minL), the IBD structure of the new SG is determined by the IBD_structure function. In this scenario, the function returns a list containing the new Grange and its IBD structure, with an action indicator set to “New”.

In essence, the `TBM_SG` function offers a flexible and adaptive methodology for SG creation and expansion, based on the characteristics of their IBD structures. It is engineered to accommodate diverse genomic patterns, rendering it a valuable tool for subsequent steps.

```{r }
## Execute the TBM_SG function with SG1 and IBDClusters as parameters (using default values for other parameters)
TBM_SG(SG1, IBDClusters)

```

## Testing the *TBM* Function:

The code then proceeds to run the TBM function. The `TBM` function extends the functionality of `TBM_SG` to generate a list of SG. It iteratively scans the chromosome using the `TBM_SG function`. The function accepts parameters such as the list of SG (which contains the first SG as recommended), IBD clusters (the only required parameter), overlap parameters, minimum length of SGs, the number of iterations N (for testing purposes), and the concerned chromosome (chr). The resultant list (SGList) contains Granges of each SG along with their IBD structures.

```{r }

## Execute the TBM function for N=50 SG with IBDClusters as a parameter
TBM(IBDClusters, N=50)

## Create a list with SG1 and execute the TBM function for N=50 SG with IBDClusters and the list as parameters
SGList=list(SG1)
TBM(IBDClusters, SGList, N=50)

```
These functions enable the dynamic creation and modification of SG based on their IBD structures, facilitating comprehensive analysis and exploration of genetic data. The TBM methodology provides flexibility and adaptability in segment generation, making it an invaluable tool for genetic researchers.

# Generated SGs and Disjointness in IBD Structures 

The `TBMD` function is designed to transform the non-disjoint clusters present in the output list of SGs (`SGList`) from the `TBM` function into a refined output comprising disjoint clusters. It sequentially examines each SG to determine if the IBD structure within the SG encompasses disjoint clusters. By default, it employs the `Disj_Clst_SG_Bigger` function for this purpose. The `Disj_Clst_SG_Bigger` function converts the IBD structure into disjoint clusters by retaining the clusters with the most individuals among the overlapped clusters.

The `TBMD` function also offers an option to build disjoint clusters by merging the overlapped clusters using the `Disj_Clst_SG` function. To optimize resources, the `Disj_Clst_SG` function applies principles of graph theory to convert the IBD structure into disjoint clusters. It employs graph theory concepts to represent the sharing between individuals, where an edge exists between individuals that share the segment, and vertices represent the two haplotypes of the individual. The graph is subsequently partitioned into disjoint clusters.

The input parameters for `TBMD` encompass the SG of interest and IBD clusters. The parameter `disjM` should be set to `merge` if the `Disj_Clst_SG` function is to be used. Otherwise, the parameter `disjM` is set to `bigger` by default for the `Disj_Clst_SG_Bigger` function. The function begins by processing the input data, eliminating superfluous columns, and managing missing values. The ultimate output of the `TBMD` function is a list of SGs, each containing disjoint clusters within their IBD structures. The function yields the modified SG with disjoint clusters in its IBD structure. The function also offers an option to omit SGs that lack any clusters in their IBD structure (`remove=TRUE`).


```{r }
# Generate a list of SGs using the TBM function
SGList.chr = TBM(IBDClusters, SGList, N=50)

# Use the TBMD function to separate the clusters within each SG
TBMD(SGList.chr, IBDClusters)

```

The R code example provided executes the TBM function with the parameters IBDClusters, SGList, and N=50 to generate a list of SGs (SGList.chr). This list is then passed to the TBMD function along with IBDClusters to separate the clusters within each SG.

# Computation of Most Shared SG with Disjoint IBD Clusters

The `Smsg` function is designed to compute the `Most Shared SG` for a list of SGs, the resulting dataframe (df1 and df2) provides information about SG and the computed Smsg value.

Function Parameters:

    SGListD.chr: List of SGs with disjoint IBD structurs for a given chromosome.
    Naff: Number of affected individuals (optional). 

Output:

    A dataframe containing columns for chromosome, SG identifier, selected cluster, start and end positions, and the computed Smsg value.
    If Naff is provided, an additional column (Smsg_per_aff) is included, representing the Smsg normalized by the number of affected individuals.
    
This function is particularly beneficial for genetic data analysis, offering a quantitative measure of the most shared SG. The provision to normalize Smsg by the number of affected individuals augments the interpretability of the results in scenarios with varying numbers of affected individuals.

The R code provided first applies the `TBMD` function to the list of SG (SGList.chr) to segregate the clusters within each segment. Then, the `Smsg` function is used to compute the most shared SG for the disjoint list (SGListD.chr). This is done twice, once without specifying the number of affected individuals, and once with Naff set to 19.

```{r }
# Segregate the clusters within each SG
SGListD.chr = TBMD(SGList.chr, IBDClusters)

# Compute Smsg
df1 = Smsg(SGListD.chr)
print(df1)

# The output in case Naff is given
df2 = Smsg(SGListD.chr, 19)
print(df2)

```

The subsequent lines of code are for saving the results. It is recommended to save the list of SGs in RDS format. The second line writes the computed Smsg data to a CSV file, providing a convenient and shareable format for further analysis and visualization. RDS files can be utilized to compute the IBD measure Sall, as proposed by Whittemore and Halpern (1994), for a single SG using the `Sall_SG` function or for all SGs of a chromosome using the `Sall` function. The inputs and output of Sall are analogous to those of Smsg. An example of Sall is not presented due to the memory and time constraints required by this function.

    # Save the list of SGs in RDS format
    saveRDS(SGListD.chr, file=SGListD.chr.filename)
    
    # Write the computed Smsg data to a CSV file
    write.csv(outSmsg, paste0(outFileNameWithPath, ".csv"), row.names = FALSE)


References :

Whittemore and Halpern (1994)
