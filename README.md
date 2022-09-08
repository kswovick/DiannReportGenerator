# DiannReportGenerator

The goal of DiannReportGenerator is to make a new document that contains, for each sample, the "MaxLFQ Intensity" and "Number of Peptides Quantified" for every Protein Group.

Precursors are filtered by both Library Q-Value and Protein Group Library Q-Value at 0.01 and PEP of 0.5, and keratin IDs are removed. After filtering, MaxLFQ-based quantification (<https://doi.org/10.1074/mcp.M113.031591>) is performed. The number of precursors/peptides that were quantified within a given protein group is counted and the data is then tidied to long format where every row corresponds to a protein group (reported as the FASTA accession number) and the rows contain:

-   Gene Name

-   Protein Name

-   MaxLFQ Intensity for each sample

-   Number of Peptides for each sample

## Installation

You can install the development version of DiannReportGenerator from:

<https://github.com/kswovick/DIANN-Report-Generator>

Prior to running, the 'DiaNN R Package' developed by the Demichev lab needs to be installed. The package and instructions can be found at:

<https://github.com/vdemichev/diann-rpackage>

In short, instructions are as follows:

    install.packages("devtools")
    library(devtools)
    install_github("https://github.com/vdemichev/diann-rpackage")
    library(diann)

## Instructions

    library(DiannReportGenerator)

Put the name of the main DIA-NN report in quotes, followed by whatever you wish to call the new document in quotes.

    diann_reporter("report.tsv", "new_report.tsv")
