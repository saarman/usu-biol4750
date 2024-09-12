## Section 4.3 Worked Example

<https://bookdown.org/hhwagner1/LandGenCourse_book/WE_1.html>

The data set we will use is ralu.loci

Make sure we have the packages:

if(!require(“adegenet”)) install.packages(“adegenet”)
if(!requireNamespace(“popgraph”, quietly = TRUE)) {
install.packages(c(“RgoogleMaps”, “geosphere”, “proto”, “sampling”,
“seqinr”, “spacetime”, “spdep”), dependencies=TRUE)
remotes::install\_github(“dyerlab/popgraph”) }
if(!requireNamespace(“gstudio”, quietly = TRUE))
remotes::install\_github(“dyerlab/gstudio”)

install.packages(“here”)

Load libraries

    library(adegenet)

    ## Loading required package: ade4

    ## 
    ##    /// adegenet 2.1.10 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

    library(gstudio)

    ## Warning: replacing previous import 'dplyr::union' by 'raster::union' when
    ## loading 'gstudio'

    ## Warning: replacing previous import 'dplyr::intersect' by 'raster::intersect'
    ## when loading 'gstudio'

    ## Warning: replacing previous import 'dplyr::select' by 'raster::select' when
    ## loading 'gstudio'

    ## 
    ## Attaching package: 'gstudio'

    ## The following objects are masked from 'package:adegenet':
    ## 
    ##     alleles, ploidy

    library(LandGenCourse)
    library(tibble)
    library(here)

    ## here() starts at /uufs/chpc.utah.edu/common/home/u6036559/git/usu-biol4750

    library(vcfR)

    ## 
    ##    *****       ***   vcfR   ***       *****
    ##    This is vcfR 1.15.0 
    ##      browseVignettes('vcfR') # Documentation
    ##      citation('vcfR') # Citation
    ##    *****       *****      *****       *****

    library(pinfsc50)
    library(utils)

## Import data

straight from the package after library is loaded

    data("ralu.loci")

As a .csv file:

    if(!dir.exists(paste0(here(),"/downloads"))) dir.create(paste0(here(),"/downloads"))
    file.copy(system.file("extdata", "ralu.loci.csv", package = "LandGenCourse"),
              paste0(here(), "/downloads/ralu.loci.csv"), overwrite=FALSE)

    ## [1] FALSE

## Section 4.4 R exercise Week 1
