superSeq: Determining sufficient read depth in RNA-Seq experiments
=======

The superSeq package models the relationship between statistical power and read depth in an RNA sequencing study. Our algorithm can help predict how many additional reads, if any, should be sequenced to achieve desired statistical power.

See also [superSeq: Determining sufficient sequencing depth in RNA-Seq differential expression studies](TODO).

Installation
-------------

First install the Bioconductor dependencies:

    source("http://bioconductor.org/biocLite.R")
    biocLite(c("qvalue", "limma", "edgeR", "DESeq2", "DEXSeq", "pasilla"))

Then install the [devtools](https://github.com/hadley/devtools) package, and use it to install the [subSeq](https://github.com/StoreyLab/subSeq) and superSeq packages. 

    install.packages("devtools")
    library(devtools)
    install_github("StoreyLab/subSeq")
    install_github("StoreyLab/superSeq", build_vignettes = TRUE)

Vignette
---------------------

Once you've installed the package, you can access the vignette with

    library(superSeq)
    vignette("superSeq")

If you run into a problem or have a question about the software's usage, please open a [GitHub issue](https://github.com/StoreyLab/superSeq/issues).
