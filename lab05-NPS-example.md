# Section 5.4 Exercise

<https://bookdown.org/hhwagner1/LandGenCourse_book/r-exercise-week-2.html>

## a. Load libraries

    library(LandGenCourse)
    #library(EcoGenetics)
    library(GeNetIt)

    ## Loading required package: nlme

    library(hierfstat)
    library(adegenet)

    ## Loading required package: ade4

    ## 
    ##    /// adegenet 2.1.10 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

    ## 
    ## Attaching package: 'adegenet'

    ## The following objects are masked from 'package:hierfstat':
    ## 
    ##     Hs, read.fstat

    require(gstudio)       

    ## Loading required package: gstudio

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

    ## The following object is masked from 'package:hierfstat':
    ## 
    ##     Ho

    require(dplyr)

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     collapse

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    require(tibble) 

    ## Loading required package: tibble

    require(sf)

    ## Loading required package: sf

    ## Linking to GEOS 3.10.2, GDAL 3.4.1, PROJ 8.2.1; sf_use_s2() is TRUE

    require(popgraph)

    ## Loading required package: popgraph

    require(RgoogleMaps)

    ## Loading required package: RgoogleMaps

    ## 
    ## Thank you for using RgoogleMaps!

    ## 
    ## To acknowledge our work, please cite the package:

    ##  Markus Loecher and Karl Ropkins (2015). RgoogleMaps and loa: Unleashing R
    ##   Graphics Power on Map Tiles. Journal of Statistical Software 63(4), 1-18.

    require(geosphere)

    ## Loading required package: geosphere

    require(proto)

    ## Loading required package: proto

    require(sampling)

    ## Loading required package: sampling

    ## 
    ## Attaching package: 'sampling'

    ## The following object is masked from 'package:adegenet':
    ## 
    ##     strata

    require(seqinr)

    ## Loading required package: seqinr

    ## 
    ## Attaching package: 'seqinr'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## The following object is masked from 'package:nlme':
    ## 
    ##     gls

    require(spacetime)

    ## Loading required package: spacetime

    require(spdep)

    ## Loading required package: spdep

    ## Loading required package: spData

    ## To access larger datasets in this package, install the spDataLarge
    ## package with: `install.packages('spDataLarge',
    ## repos='https://nowosad.github.io/drat/', type='source')`

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## 
    ## Attaching package: 'spdep'

    ## The following object is masked from 'package:ade4':
    ## 
    ##     mstree

    require(here)

    ## Loading required package: here

    ## here() starts at /uufs/chpc.utah.edu/common/home/u6036559/git/usu-biol4750

1.  Import to adegenet object with gstudio

<!-- -->

    library(gstudio)
    library(adegenet)

    # 1. CSV file "./downloads/pulsatilla_genotypes.csv" --> data frame 
    # with "gstudio" function read_population()
    g.Flr <- read_population("./downloads/pulsatilla_genotypes.csv",
                           type = "column",locus.columns = c(6:19))

    # 2. Select only adults with base R indexing of data frame 
    # rows where OffID==0, all columns
    g.Flr <- g.Flr[g.Flr$OffID==0,]

    # 3. Nothing to do here

    # 4. Create genind object with "adegenet" function df2genind() 
    # using NA.char = ""
    g.Flr.genind <- df2genind(X=g.Flr[,c(6:12)], sep=":", ncode=NULL, ind.names=g.Flr$ID, loc.names=NULL, pop=g.Flr$Population, NA.char="", ploidy=2, type="codom", strata=NULL, hierarchy=NULL)

    # 5. Check genind object
    g.Flr.genind

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 221 individuals; 7 loci; 105 alleles; size: 129.8 Kb
    ## 
    ##  // Basic content
    ##    @tab:  221 x 105 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 8-25)
    ##    @loc.fac: locus factor for the 105 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = g.Flr[, c(6:12)], sep = ":", ncode = NULL, ind.names = g.Flr$ID, 
    ##     loc.names = NULL, pop = g.Flr$Population, NA.char = "", ploidy = 2, 
    ##     type = "codom", strata = NULL, hierarchy = NULL)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 14-56)

    summary(g.Flr.genind)

    ## 
    ## // Number of individuals: 221
    ## // Group sizes: 21 56 21 22 14 42 45
    ## // Number of alleles per locus: 18 8 25 8 19 14 13
    ## // Number of alleles per group: 63 68 54 50 51 73 53
    ## // Percentage of missing data: 0.9 %
    ## // Observed heterozygosity: 0.74 0.54 0.89 0.71 0.74 0.68 0.74
    ## // Expected heterozygosity: 0.83 0.57 0.89 0.74 0.81 0.76 0.83

Without G-Studio, Drop offspring Using adegenet and base R

    library(adegenet)

    # 1. CSV file "./downloads/pulsatilla_genotypes.csv" --> data frame 
    # with base R function read.csv()
    Flr <- read.csv("./downloads/pulsatilla_genotypes.csv", header=TRUE)

    # 2. Select only adults with base R indexing of data frame 
    # rows where OffID==0, all columns
    Flr <- Flr[Flr$OffID==0,]

    # 3. Combine columns with base R function paste()
    Flr <- data.frame(Flr[,1:5],loc1 = paste(Flr[,6],  Flr[,7], sep=":"), 
                                loc2 = paste(Flr[,8],  Flr[,9], sep=":"), 
                                loc3 = paste(Flr[,10], Flr[,11], sep=":"), 
                                loc4 = paste(Flr[,12], Flr[,13], sep=":"), 
                                loc5 = paste(Flr[,14], Flr[,15], sep=":"),  
                                loc6 = paste(Flr[,16], Flr[,17], sep=":"),  
                                loc7 = paste(Flr[,18], Flr[,19], sep=":"))
    # 4. Create genind object with "adegenet" function df2genind() 
    # using NA.char = "NA"
    Flr.genind <- df2genind(X=Flr[,c(6:12)], sep=":", ncode=NULL, ind.names= Flr$ID, loc.names=names(Flr[,c(6:12)]), pop=Flr$Population, NA.char="NA", ploidy=2, type="codom", strata=NULL, hierarchy=NULL)

    # 5. Check genind object
    Flr.genind

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 221 individuals; 7 loci; 105 alleles; size: 130.8 Kb
    ## 
    ##  // Basic content
    ##    @tab:  221 x 105 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 8-25)
    ##    @loc.fac: locus factor for the 105 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = Flr[, c(6:12)], sep = ":", ncode = NULL, ind.names = Flr$ID, 
    ##     loc.names = names(Flr[, c(6:12)]), pop = Flr$Population, 
    ##     NA.char = "NA", ploidy = 2, type = "codom", strata = NULL, 
    ##     hierarchy = NULL)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 14-56)

    Flr.genind@pop

    ##   [1] A21  A21  A21  A21  A21  A21  A21  A21  A21  A21  A21  A21  A21  A21  A21 
    ##  [16] A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25 
    ##  [31] A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25 
    ##  [46] A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25  A25 
    ##  [61] A25  A25  A25  A26  A26  A26  A26  A26  A26  A26  A26  A26  A26  A26  A26 
    ##  [76] A26  A26  A45  A45  A45  A45  A45  A45  A45  A45  A45  A45  A45  A45  A45 
    ##  [91] A45  A45  A45  A41  A41  A41  A41  A41  A03  A03  A03  A03  A03  A03  A03 
    ## [106] A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03 
    ## [121] A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03  A03 
    ## [136] G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a
    ## [151] G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a G05a
    ## [166] G05a G05a G05a G05a G05a G05a G05a A25  A25  A25  A25  A25  A25  A25  A25 
    ## [181] A26  A26  A26  A26  A26  A26  A26  A45  A45  A45  A45  A45  A45  A41  A41 
    ## [196] A41  A41  A41  A41  A41  A41  A41  A21  A21  A21  A21  A21  A21  G05a G05a
    ## [211] G05a G05a G05a G05a G05a A03  A03  A03  A03  A03  G05a
    ## Levels: A21 A25 A26 A45 A41 A03 G05a

    summary(Flr.genind)

    ## 
    ## // Number of individuals: 221
    ## // Group sizes: 21 56 21 22 14 42 45
    ## // Number of alleles per locus: 18 8 25 8 19 14 13
    ## // Number of alleles per group: 63 68 54 50 51 73 53
    ## // Percentage of missing data: 0.9 %
    ## // Observed heterozygosity: 0.74 0.54 0.89 0.71 0.74 0.68 0.74
    ## // Expected heterozygosity: 0.83 0.57 0.89 0.74 0.81 0.76 0.83

1.  PLot locations of individuals from site A25 From dataframe

<!-- -->

    # Select data for Population "A25"
    Sites <- Flr[Flr$Population == "A25", c("X", "Y")]

    # Check if the data frame has valid coordinates
    head(Sites)

    ##          X       Y
    ## 16 4422658 5425371
    ## 17 4422659 5425372
    ## 18 4422659 5425371
    ## 19 4422658 5425370
    ## 20 4422659 5425371
    ## 21 4422659 5425371

    # Convert to an 'sf' object
    Pulsatilla <- st_as_sf(Sites, coords = c("X", "Y"), crs = 31468) 

    # Confirm the CRS was set correctly
    st_crs(Pulsatilla)

    ## Coordinate Reference System:
    ##   User input: EPSG:31468 
    ##   wkt:
    ## PROJCRS["DHDN / 3-degree Gauss-Kruger zone 4",
    ##     BASEGEOGCRS["DHDN",
    ##         DATUM["Deutsches Hauptdreiecksnetz",
    ##             ELLIPSOID["Bessel 1841",6377397.155,299.1528128,
    ##                 LENGTHUNIT["metre",1]]],
    ##         PRIMEM["Greenwich",0,
    ##             ANGLEUNIT["degree",0.0174532925199433]],
    ##         ID["EPSG",4314]],
    ##     CONVERSION["3-degree Gauss-Kruger zone 4",
    ##         METHOD["Transverse Mercator",
    ##             ID["EPSG",9807]],
    ##         PARAMETER["Latitude of natural origin",0,
    ##             ANGLEUNIT["degree",0.0174532925199433],
    ##             ID["EPSG",8801]],
    ##         PARAMETER["Longitude of natural origin",12,
    ##             ANGLEUNIT["degree",0.0174532925199433],
    ##             ID["EPSG",8802]],
    ##         PARAMETER["Scale factor at natural origin",1,
    ##             SCALEUNIT["unity",1],
    ##             ID["EPSG",8805]],
    ##         PARAMETER["False easting",4500000,
    ##             LENGTHUNIT["metre",1],
    ##             ID["EPSG",8806]],
    ##         PARAMETER["False northing",0,
    ##             LENGTHUNIT["metre",1],
    ##             ID["EPSG",8807]]],
    ##     CS[Cartesian,2],
    ##         AXIS["northing (X)",north,
    ##             ORDER[1],
    ##             LENGTHUNIT["metre",1]],
    ##         AXIS["easting (Y)",east,
    ##             ORDER[2],
    ##             LENGTHUNIT["metre",1]],
    ##     USAGE[
    ##         SCOPE["Cadastre, engineering survey, topographic mapping."],
    ##         AREA["Germany - former West Germany onshore between 10°30'E and 13°30'E - states of Bayern, Berlin, Niedersachsen, Schleswig-Holstein."],
    ##         BBOX[47.39,10.5,54.59,13.51]],
    ##     ID["EPSG",31468]]

    # Transform to WGS84 (lat/long, EPSG:4326)
    Pulsatilla_latlon <- st_transform(Pulsatilla, crs = 4326)

# d. Add geographic info to genind object

Genind object @other can hold a list such as spatial coordinates

    # Convert data frame of xy coordinates from Flr into a matrix of xy coordinates
    xy_matrix <- as.matrix(Flr[, c("X", "Y")])

    # Add the matrix to the @other slot with the name "xy"
    Flr.genind@other <- list(xy = xy_matrix)

    # Check the class of the @other slot
    class(Flr.genind@other)  # Should return "list"

    ## [1] "list"

    # Inspect the genind object and its summary
    print(Flr.genind)

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 221 individuals; 7 loci; 105 alleles; size: 148.9 Kb
    ## 
    ##  // Basic content
    ##    @tab:  221 x 105 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 8-25)
    ##    @loc.fac: locus factor for the 105 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = Flr[, c(6:12)], sep = ":", ncode = NULL, ind.names = Flr$ID, 
    ##     loc.names = names(Flr[, c(6:12)]), pop = Flr$Population, 
    ##     NA.char = "NA", ploidy = 2, type = "codom", strata = NULL, 
    ##     hierarchy = NULL)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 14-56)
    ##    @other: a list containing: xy

    summary(Flr.genind)

    ## 
    ## // Number of individuals: 221
    ## // Group sizes: 21 56 21 22 14 42 45
    ## // Number of alleles per locus: 18 8 25 8 19 14 13
    ## // Number of alleles per group: 63 68 54 50 51 73 53
    ## // Percentage of missing data: 0.9 %
    ## // Observed heterozygosity: 0.74 0.54 0.89 0.71 0.74 0.68 0.74
    ## // Expected heterozygosity: 0.83 0.57 0.89 0.74 0.81 0.76 0.83

# e. Calculate genetic and geographic euclidean distance (for Site A25)

Genetic with adegenet::propShared

    #Check genind object
    Flr.genind

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 221 individuals; 7 loci; 105 alleles; size: 148.9 Kb
    ## 
    ##  // Basic content
    ##    @tab:  221 x 105 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 8-25)
    ##    @loc.fac: locus factor for the 105 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = Flr[, c(6:12)], sep = ":", ncode = NULL, ind.names = Flr$ID, 
    ##     loc.names = names(Flr[, c(6:12)]), pop = Flr$Population, 
    ##     NA.char = "NA", ploidy = 2, type = "codom", strata = NULL, 
    ##     hierarchy = NULL)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 14-56)
    ##    @other: a list containing: xy

    summary(Flr.genind)

    ## 
    ## // Number of individuals: 221
    ## // Group sizes: 21 56 21 22 14 42 45
    ## // Number of alleles per locus: 18 8 25 8 19 14 13
    ## // Number of alleles per group: 63 68 54 50 51 73 53
    ## // Percentage of missing data: 0.9 %
    ## // Observed heterozygosity: 0.74 0.54 0.89 0.71 0.74 0.68 0.74
    ## // Expected heterozygosity: 0.83 0.57 0.89 0.74 0.81 0.76 0.83

    # Subset only site A25
    Flr.A25 <- Flr.genind[Flr.genind@pop == "A25"]

    # Calculate genetic distance
    ind.dist <- adegenet::propShared(Flr.A25)
    ind.dist <- as.dist(1-ind.dist)
    class(ind.dist)

    ## [1] "dist"

Geographic with dist()

    # We have already changed coordinates and subset A25
    # Extract coordinates from the sf object as a matrix
    coords_wgs84 <- st_coordinates(Pulsatilla_latlon)

    # Print the extracted WGS84 coordinates to verify
    print(coords_wgs84)

    ##              X        Y
    ##  [1,] 10.94236 48.96085
    ##  [2,] 10.94239 48.96086
    ##  [3,] 10.94238 48.96086
    ##  [4,] 10.94237 48.96085
    ##  [5,] 10.94238 48.96086
    ##  [6,] 10.94239 48.96086
    ##  [7,] 10.94238 48.96086
    ##  [8,] 10.94237 48.96085
    ##  [9,] 10.94239 48.96087
    ## [10,] 10.94239 48.96084
    ## [11,] 10.94238 48.96084
    ## [12,] 10.94237 48.96083
    ## [13,] 10.94240 48.96082
    ## [14,] 10.94238 48.96080
    ## [15,] 10.94238 48.96081
    ## [16,] 10.94238 48.96081
    ## [17,] 10.94239 48.96080
    ## [18,] 10.94241 48.96079
    ## [19,] 10.94238 48.96081
    ## [20,] 10.94239 48.96079
    ## [21,] 10.94237 48.96080
    ## [22,] 10.94236 48.96080
    ## [23,] 10.94241 48.96077
    ## [24,] 10.94238 48.96079
    ## [25,] 10.94239 48.96079
    ## [26,] 10.94239 48.96078
    ## [27,] 10.94238 48.96077
    ## [28,] 10.94234 48.96076
    ## [29,] 10.94231 48.96075
    ## [30,] 10.94233 48.96075
    ## [31,] 10.94235 48.96076
    ## [32,] 10.94233 48.96077
    ## [33,] 10.94231 48.96070
    ## [34,] 10.94232 48.96070
    ## [35,] 10.94232 48.96061
    ## [36,] 10.94235 48.96063
    ## [37,] 10.94233 48.96062
    ## [38,] 10.94232 48.96062
    ## [39,] 10.94232 48.96064
    ## [40,] 10.94232 48.96064
    ## [41,] 10.94233 48.96063
    ## [42,] 10.94233 48.96063
    ## [43,] 10.94233 48.96063
    ## [44,] 10.94235 48.96064
    ## [45,] 10.94234 48.96064
    ## [46,] 10.94250 48.96070
    ## [47,] 10.94253 48.96072
    ## [48,] 10.94253 48.96072
    ## [49,] 10.94235 48.96086
    ## [50,] 10.94245 48.96080
    ## [51,] 10.94243 48.96080
    ## [52,] 10.94238 48.96085
    ## [53,] 10.94239 48.96085
    ## [54,] 10.94241 48.96083
    ## [55,] 10.94238 48.96078
    ## [56,] 10.94241 48.96085

    # Calculate the Euclidean distance using dist()
    euclidean_distances <- dist(coords_wgs84, method = "euclidean")

    # Print the Euclidean distance matrix
    print(euclidean_distances)

    ##               1            2            3            4            5
    ## 2  2.228977e-05                                                    
    ## 3  1.923686e-05 4.634166e-06                                       
    ## 4  6.300033e-06 1.979748e-05 1.587047e-05                          
    ## 5  1.922449e-05 5.807367e-06 1.351288e-06 1.545900e-05             
    ## 6  2.510122e-05 4.476402e-06 5.882784e-06 2.165374e-05 6.210366e-06
    ## 7  1.947767e-05 7.354475e-06 3.084638e-06 1.522315e-05 1.735537e-06
    ## 8  9.439574e-06 1.911386e-05 1.481852e-05 3.140692e-06 1.416699e-05
    ## 9  2.991668e-05 7.861282e-06 1.215316e-05 2.765858e-05 1.304162e-05
    ## 10 2.277153e-05 1.748207e-05 1.340338e-05 1.670609e-05 1.205229e-05
    ## 11 2.362780e-05 2.710495e-05 2.259257e-05 1.755577e-05 2.132125e-05
    ## 12 2.527583e-05 3.762699e-05 3.300297e-05 2.132719e-05 3.194465e-05
    ## 13 4.379720e-05 4.005196e-05 3.662231e-05 3.754365e-05 3.529671e-05
    ## 14 5.133393e-05 5.835646e-05 5.393444e-05 4.652001e-05 5.263378e-05
    ## 15 4.495374e-05 4.931212e-05 4.503152e-05 3.953313e-05 4.370153e-05
    ## 16 4.993929e-05 5.504723e-05 5.074683e-05 4.472129e-05 4.942026e-05
    ## 17 5.890550e-05 5.992038e-05 5.606701e-05 5.315132e-05 5.471612e-05
    ## 18 7.975462e-05 7.719937e-05 7.395145e-05 7.373090e-05 7.263601e-05
    ## 19 5.000466e-05 5.516987e-05 5.086564e-05 4.479810e-05 4.953968e-05
    ## 20 6.358186e-05 6.880090e-05 6.456305e-05 5.852670e-05 6.322866e-05
    ## 21 5.779112e-05 6.937987e-05 6.477550e-05 5.425869e-05 6.357374e-05
    ## 22 5.272705e-05 6.609410e-05 6.146043e-05 4.974139e-05 6.032889e-05
    ## 23 9.205171e-05 9.268110e-05 8.904030e-05 8.641622e-05 8.769526e-05
    ## 24 6.345088e-05 7.069919e-05 6.630927e-05 5.881536e-05 6.500099e-05
    ## 25 7.194949e-05 7.691207e-05 7.272897e-05 6.690757e-05 7.138914e-05
    ## 26 7.801185e-05 8.262707e-05 7.849629e-05 7.294185e-05 7.715241e-05
    ## 27 8.035659e-05 8.756099e-05 8.321945e-05 7.582839e-05 8.190145e-05
    ## 28 9.360660e-05 1.092390e-04 1.046115e-04 9.173190e-05 1.035327e-04
    ## 29 1.178172e-04 1.358390e-04 1.312923e-04 1.170912e-04 1.303404e-04
    ## 30 1.042354e-04 1.200887e-04 1.154627e-04 1.024812e-04 1.143886e-04
    ## 31 9.608249e-05 1.100113e-04 1.053817e-04 9.355751e-05 1.042265e-04
    ## 32 9.353299e-05 1.106091e-04 1.060222e-04 9.229639e-05 1.050223e-04
    ## 33 1.633966e-04 1.794929e-04 1.748648e-04 1.618349e-04 1.737823e-04
    ## 34 1.588487e-04 1.739156e-04 1.692816e-04 1.568704e-04 1.681511e-04
    ## 35 2.493735e-04 2.624470e-04 2.578582e-04 2.467617e-04 2.566399e-04
    ## 36 2.285705e-04 2.394862e-04 2.349961e-04 2.252857e-04 2.337204e-04
    ## 37 2.347088e-04 2.468028e-04 2.422523e-04 2.317773e-04 2.410058e-04
    ## 38 2.342328e-04 2.473729e-04 2.427804e-04 2.316340e-04 2.415655e-04
    ## 39 2.205819e-04 2.338679e-04 2.292695e-04 2.180226e-04 2.280607e-04
    ## 40 2.186471e-04 2.316381e-04 2.270492e-04 2.159890e-04 2.258310e-04
    ## 41 2.273238e-04 2.399628e-04 2.353880e-04 2.245581e-04 2.341580e-04
    ## 42 2.223938e-04 2.350900e-04 2.305124e-04 2.196430e-04 2.292846e-04
    ## 43 2.260161e-04 2.382472e-04 2.336893e-04 2.231212e-04 2.324475e-04
    ## 44 2.132664e-04 2.238215e-04 2.193496e-04 2.098682e-04 2.180671e-04
    ## 45 2.169744e-04 2.283546e-04 2.238365e-04 2.138157e-04 2.225728e-04
    ## 46 2.009674e-04 1.949436e-04 1.925562e-04 1.948147e-04 1.913418e-04
    ## 47 2.172996e-04 2.079031e-04 2.062146e-04 2.110107e-04 2.051170e-04
    ## 48 2.117678e-04 2.022267e-04 2.005655e-04 2.054760e-04 1.994732e-04
    ## 49 1.134716e-05 3.208976e-05 2.977788e-05 1.734842e-05 2.998493e-05
    ## 50 9.879039e-05 8.808968e-05 8.648808e-05 9.250326e-05 8.542075e-05
    ## 51 8.878308e-05 8.011840e-05 7.801800e-05 8.248478e-05 7.686154e-05
    ## 52 1.787987e-05 1.238559e-05 7.843749e-06 1.247902e-05 6.582884e-06
    ## 53 2.206941e-05 1.005754e-05 6.571021e-06 1.711843e-05 5.316563e-06
    ## 54 4.627957e-05 3.692360e-05 3.452570e-05 4.003693e-05 3.334339e-05
    ## 55 7.379809e-05 8.188837e-05 7.747103e-05 6.943608e-05 7.617010e-05
    ## 56 4.568302e-05 2.935771e-05 2.888854e-05 4.015700e-05 2.813694e-05
    ##               6            7            8            9           10
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7  6.889024e-06                                                    
    ## 8  2.035136e-05 1.360905e-05                                       
    ## 9  7.731828e-06 1.419886e-05 2.691281e-05                          
    ## 10 1.542617e-05 1.032568e-05 1.372341e-05 2.313776e-05             
    ## 11 2.579653e-05 1.975868e-05 1.464167e-05 3.352128e-05 1.072606e-05
    ## 12 3.742227e-05 3.070988e-05 1.991348e-05 4.490864e-05 2.367035e-05
    ## 13 3.675963e-05 3.357498e-05 3.444919e-05 4.363029e-05 2.366969e-05
    ## 14 5.651632e-05 5.100937e-05 4.433511e-05 6.419552e-05 4.110375e-05
    ## 15 4.715918e-05 4.202502e-05 3.699721e-05 5.476207e-05 3.186771e-05
    ## 16 5.291297e-05 4.774986e-05 4.229703e-05 6.051442e-05 3.761380e-05
    ## 17 5.701753e-05 5.298294e-05 5.038032e-05 6.417512e-05 4.267786e-05
    ## 18 7.363234e-05 7.092237e-05 7.078069e-05 7.996501e-05 6.103481e-05
    ## 19 5.304370e-05 4.787043e-05 4.238065e-05 6.064776e-05 3.773975e-05
    ## 20 6.649227e-05 6.154252e-05 5.616966e-05 7.400860e-05 5.132601e-05
    ## 21 6.829219e-05 6.210599e-05 5.278329e-05 7.601548e-05 5.305580e-05
    ## 22 6.543471e-05 5.896774e-05 4.858792e-05 7.308961e-05 5.057397e-05
    ## 23 8.946450e-05 8.596012e-05 8.369341e-05 9.624695e-05 7.573512e-05
    ## 24 6.873588e-05 6.336242e-05 5.670186e-05 7.637708e-05 5.336259e-05
    ## 25 7.448013e-05 6.969174e-05 6.454528e-05 8.192701e-05 5.943020e-05
    ## 26 8.009102e-05 7.544606e-05 7.055597e-05 8.747282e-05 6.515501e-05
    ## 27 8.544449e-05 8.024514e-05 7.374617e-05 9.302602e-05 7.014805e-05
    ## 28 1.088173e-04 1.022414e-04 9.104761e-05 1.164199e-04 9.407938e-05
    ## 29 1.360201e-04 1.292264e-04 1.169585e-04 1.433794e-04 1.220203e-04
    ## 30 1.196855e-04 1.131028e-04 1.018431e-04 1.272841e-04 1.049488e-04
    ## 31 1.091626e-04 1.028237e-04 9.253182e-05 1.168634e-04 9.403058e-05
    ## 32 1.105773e-04 1.038435e-04 9.193867e-05 1.180344e-04 9.632575e-05
    ## 33 1.790266e-04 1.724812e-04 1.612449e-04 1.866498e-04 1.641672e-04
    ## 34 1.731832e-04 1.667809e-04 1.560703e-04 1.808689e-04 1.580956e-04
    ## 35 2.611333e-04 2.551390e-04 2.456060e-04 2.688625e-04 2.457389e-04
    ## 36 2.377175e-04 2.321317e-04 2.237874e-04 2.453820e-04 2.223084e-04
    ## 37 2.452798e-04 2.394626e-04 2.304603e-04 2.529904e-04 2.298551e-04
    ## 38 2.460860e-04 2.400700e-04 2.304884e-04 2.538164e-04 2.306991e-04
    ## 39 2.326265e-04 2.265744e-04 2.169007e-04 2.403581e-04 2.172535e-04
    ## 40 2.303302e-04 2.243306e-04 2.148169e-04 2.380599e-04 2.149387e-04
    ## 41 2.385671e-04 2.326397e-04 2.233283e-04 2.462914e-04 2.231568e-04
    ## 42 2.337116e-04 2.277697e-04 2.184222e-04 2.414373e-04 2.183045e-04
    ## 43 2.367617e-04 2.309112e-04 2.218251e-04 2.444772e-04 2.213398e-04
    ## 44 2.219931e-04 2.164677e-04 2.083146e-04 2.296429e-04 2.065963e-04
    ## 45 2.266935e-04 2.210030e-04 2.123855e-04 2.343820e-04 2.112703e-04
    ## 46 1.908009e-04 1.897309e-04 1.917758e-04 1.956148e-04 1.806988e-04
    ## 47 2.035331e-04 2.036382e-04 2.078839e-04 2.072450e-04 1.956021e-04
    ## 48 1.978509e-04 1.980007e-04 2.023468e-04 2.015257e-04 1.900139e-04
    ## 49 3.548724e-05 3.046411e-05 2.042784e-05 3.917446e-05 3.403200e-05
    ## 50 8.371725e-05 8.398336e-05 8.936434e-05 8.758912e-05 7.645180e-05
    ## 51 7.589603e-05 7.532015e-05 7.934948e-05 8.058404e-05 6.700269e-05
    ## 52 1.174522e-05 5.105805e-06 1.010345e-05 1.925786e-05 6.445302e-06
    ## 53 7.872767e-06 3.733758e-06 1.495067e-05 1.560433e-05 7.600988e-06
    ## 54 3.291025e-05 3.178122e-05 3.691764e-05 3.860454e-05 2.376715e-05
    ## 55 7.997838e-05 7.454347e-05 6.745706e-05 8.762856e-05 6.459089e-05
    ## 56 2.488176e-05 2.712776e-05 3.743797e-05 2.801371e-05 2.411001e-05
    ##              11           12           13           14           15
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12 1.361466e-05                                                    
    ## 13 2.053007e-05 2.889268e-05                                       
    ## 14 3.142206e-05 2.648400e-05 2.674694e-05                          
    ## 15 2.301331e-05 2.201844e-05 1.703275e-05 1.008288e-05             
    ## 16 2.858688e-05 2.601524e-05 2.157815e-05 5.350903e-06 5.753830e-06
    ## 17 3.581168e-05 3.688128e-05 2.087311e-05 1.631442e-05 1.488493e-05
    ## 18 5.617537e-05 5.894198e-05 3.737897e-05 3.710839e-05 3.704343e-05
    ## 19 2.869125e-05 2.604427e-05 2.174577e-05 5.185260e-06 5.885778e-06
    ## 20 4.252094e-05 3.899791e-05 3.263002e-05 1.259314e-05 1.955683e-05
    ## 21 4.249568e-05 3.293790e-05 4.235744e-05 1.597785e-05 2.533388e-05
    ## 22 3.985247e-05 2.871909e-05 4.321991e-05 1.883869e-05 2.647245e-05
    ## 23 6.916020e-05 6.878786e-05 5.270494e-05 4.325313e-05 4.720743e-05
    ## 24 4.384042e-05 3.831696e-05 3.671865e-05 1.243707e-05 2.167175e-05
    ## 25 5.082676e-05 4.726984e-05 3.969261e-05 2.079056e-05 2.782207e-05
    ## 26 5.674415e-05 5.334221e-05 4.475237e-05 2.686050e-05 3.373139e-05
    ## 27 6.083073e-05 5.512842e-05 5.178961e-05 2.948004e-05 3.828953e-05
    ## 28 8.335429e-05 7.161460e-05 8.447799e-05 5.776269e-05 6.752485e-05
    ## 29 1.113838e-04 9.864554e-05 1.148911e-04 8.841107e-05 9.785838e-05
    ## 30 9.422401e-05 8.246201e-05 9.501751e-05 6.827261e-05 7.812457e-05
    ## 31 8.339863e-05 7.270141e-05 8.194077e-05 5.533012e-05 6.537687e-05
    ## 32 8.564492e-05 7.315504e-05 8.893604e-05 6.256051e-05 7.190507e-05
    ## 33 1.534564e-04 1.418691e-04 1.525941e-04 1.260482e-04 1.361012e-04
    ## 34 1.474467e-04 1.364121e-04 1.452195e-04 1.189735e-04 1.290564e-04
    ## 35 2.354174e-04 2.257003e-04 2.294777e-04 2.049636e-04 2.148109e-04
    ## 36 2.124035e-04 2.039694e-04 2.041716e-04 1.812047e-04 1.906751e-04
    ## 37 2.197046e-04 2.105560e-04 2.127337e-04 1.888524e-04 1.985541e-04
    ## 38 2.203567e-04 2.105851e-04 2.146144e-04 1.899726e-04 1.998432e-04
    ## 39 2.068774e-04 1.970040e-04 2.014342e-04 1.766116e-04 1.865157e-04
    ## 40 2.046103e-04 1.949087e-04 1.988283e-04 1.741921e-04 1.840578e-04
    ## 41 2.128971e-04 2.034149e-04 2.066175e-04 1.822854e-04 1.920881e-04
    ## 42 2.080307e-04 1.985090e-04 2.018621e-04 1.774572e-04 1.872745e-04
    ## 43 2.111556e-04 2.019166e-04 2.044120e-04 1.803721e-04 1.901081e-04
    ## 44 1.967603e-04 1.885422e-04 1.882607e-04 1.654982e-04 1.749057e-04
    ## 45 2.012521e-04 1.925236e-04 1.936221e-04 1.701878e-04 1.797729e-04
    ## 46 1.773441e-04 1.803375e-04 1.574542e-04 1.562681e-04 1.583251e-04
    ## 47 1.940676e-04 1.992421e-04 1.735822e-04 1.773879e-04 1.776458e-04
    ## 48 1.885756e-04 1.938841e-04 1.680773e-04 1.722303e-04 1.723428e-04
    ## 49 3.353419e-05 3.142978e-05 5.404169e-05 5.786420e-05 5.289157e-05
    ## 50 7.648007e-05 8.436355e-05 5.617455e-05 6.918273e-05 6.522365e-05
    ## 51 6.585257e-05 7.278239e-05 4.533351e-05 5.652658e-05 5.296801e-05
    ## 52 1.474882e-05 2.568458e-05 3.010361e-05 4.610970e-05 3.730958e-05
    ## 53 1.794438e-05 2.995982e-05 3.013979e-05 4.870173e-05 3.944846e-05
    ## 54 2.551836e-05 3.686228e-05 1.097786e-05 3.767936e-05 2.800946e-05
    ## 55 5.494365e-05 4.852798e-05 4.760792e-05 2.353660e-05 3.293498e-05
    ## 56 3.149395e-05 4.490828e-05 2.561799e-05 5.205335e-05 4.204360e-05
    ##              16           17           18           19           20
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17 1.252744e-05                                                    
    ## 18 3.440922e-05 2.220514e-05                                       
    ## 19 1.685925e-07 1.260327e-05 3.445594e-05                          
    ## 20 1.393972e-05 1.378528e-05 2.808373e-05 1.384101e-05             
    ## 21 2.132602e-05 3.073452e-05 4.827086e-05 2.116124e-05 2.018727e-05
    ## 22 2.381926e-05 3.490672e-05 5.403574e-05 2.367269e-05 2.617737e-05
    ## 23 4.277400e-05 3.335740e-05 1.856063e-05 4.274799e-05 3.102677e-05
    ## 24 1.596461e-05 1.957368e-05 3.410193e-05 1.581958e-05 6.375162e-06
    ## 25 2.228094e-05 1.931541e-05 2.613684e-05 2.218840e-05 8.380884e-06
    ## 26 2.826431e-05 2.397888e-05 2.563523e-05 2.817765e-05 1.443058e-05
    ## 27 3.253762e-05 3.186170e-05 3.616773e-05 3.240985e-05 1.923669e-05
    ## 28 6.303828e-05 6.975659e-05 8.030268e-05 6.286974e-05 5.626692e-05
    ## 29 9.374852e-05 1.011195e-04 1.115753e-04 9.358105e-05 8.772618e-05
    ## 30 7.350402e-05 7.966078e-05 8.877929e-05 7.333546e-05 6.600726e-05
    ## 31 6.037200e-05 6.538714e-05 7.344369e-05 6.020511e-05 5.161856e-05
    ## 32 6.790837e-05 7.580717e-05 8.800902e-05 6.774172e-05 6.272325e-05
    ## 33 1.310524e-04 1.351717e-04 1.390402e-04 1.308862e-04 1.214338e-04
    ## 34 1.238163e-04 1.270781e-04 1.295899e-04 1.236528e-04 1.134583e-04
    ## 35 2.091801e-04 2.094245e-04 2.053919e-04 2.090299e-04 1.968490e-04
    ## 36 1.849378e-04 1.835675e-04 1.769654e-04 1.847991e-04 1.718413e-04
    ## 37 1.928642e-04 1.924070e-04 1.873041e-04 1.927187e-04 1.801703e-04
    ## 38 1.942265e-04 1.946487e-04 1.910675e-04 1.940754e-04 1.819848e-04
    ## 39 1.809214e-04 1.815964e-04 1.786033e-04 1.807690e-04 1.688145e-04
    ## 40 1.784391e-04 1.788802e-04 1.755064e-04 1.782880e-04 1.661990e-04
    ## 41 1.864364e-04 1.864958e-04 1.823618e-04 1.862876e-04 1.739962e-04
    ## 42 1.816299e-04 1.817832e-04 1.778674e-04 1.814806e-04 1.692361e-04
    ## 43 1.844299e-04 1.841588e-04 1.794701e-04 1.842833e-04 1.718214e-04
    ## 44 1.691620e-04 1.676221e-04 1.609149e-04 1.690249e-04 1.559903e-04
    ## 45 1.740544e-04 1.731574e-04 1.675020e-04 1.739123e-04 1.611488e-04
    ## 46 1.548984e-04 1.434572e-04 1.215565e-04 1.549082e-04 1.444040e-04
    ## 47 1.751158e-04 1.629137e-04 1.407376e-04 1.751560e-04 1.664778e-04
    ## 48 1.698859e-04 1.576428e-04 1.354885e-04 1.699285e-04 1.614245e-04
    ## 49 5.734194e-05 6.743346e-05 8.892661e-05 5.738156e-05 7.041433e-05
    ## 50 6.500299e-05 5.286957e-05 3.610832e-05 6.511002e-05 6.296919e-05
    ## 51 5.245854e-05 4.021478e-05 2.383208e-05 5.256031e-05 5.023003e-05
    ## 52 4.300242e-05 4.880301e-05 6.745447e-05 4.311802e-05 5.686161e-05
    ## 53 4.519813e-05 4.986454e-05 6.742547e-05 4.532541e-05 5.887138e-05
    ## 54 3.245619e-05 2.978482e-05 4.165946e-05 3.262455e-05 4.264121e-05
    ## 55 2.721984e-05 2.889138e-05 3.763713e-05 2.707655e-05 1.524817e-05
    ## 56 4.704496e-05 4.539480e-05 5.591954e-05 4.721037e-05 5.806252e-05
    ##              21           22           23           24           25
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17                                                                 
    ## 18                                                                 
    ## 19                                                                 
    ## 20                                                                 
    ## 21                                                                 
    ## 22 7.247707e-06                                                    
    ## 23 4.883601e-05 5.578037e-05                                       
    ## 24 1.445031e-05 2.102631e-05 3.487134e-05                          
    ## 25 2.460931e-05 3.147213e-05 2.431858e-05 1.058510e-05             
    ## 26 2.936840e-05 3.646100e-05 1.969953e-05 1.612595e-05 6.074456e-06
    ## 27 2.595867e-05 3.315690e-05 2.771667e-05 1.705405e-05 1.267660e-05
    ## 28 4.227844e-05 4.351506e-05 7.132425e-05 5.019024e-05 5.455946e-05
    ## 29 7.255087e-05 7.240161e-05 1.014230e-04 8.158878e-05 8.601126e-05
    ## 30 5.297386e-05 5.438098e-05 7.828636e-05 6.010480e-05 6.355754e-05
    ## 31 4.103952e-05 4.398837e-05 6.287924e-05 4.601088e-05 4.852586e-05
    ## 32 4.663827e-05 4.641975e-05 8.008030e-05 5.644744e-05 6.197053e-05
    ## 33 1.115499e-04 1.136172e-04 1.241028e-04 1.162600e-04 1.169237e-04
    ## 34 1.051085e-04 1.078225e-04 1.140630e-04 1.085968e-04 1.084237e-04
    ## 35 1.931072e-04 1.970366e-04 1.872817e-04 1.932480e-04 1.901134e-04
    ## 36 1.710316e-04 1.758493e-04 1.585039e-04 1.690179e-04 1.644791e-04
    ## 37 1.777225e-04 1.820709e-04 1.690257e-04 1.768996e-04 1.731494e-04
    ## 38 1.780194e-04 1.819139e-04 1.730733e-04 1.783085e-04 1.753334e-04
    ## 39 1.644968e-04 1.683175e-04 1.607723e-04 1.650294e-04 1.622848e-04
    ## 40 1.622996e-04 1.662563e-04 1.575900e-04 1.625198e-04 1.595648e-04
    ## 41 1.707004e-04 1.748170e-04 1.642657e-04 1.704878e-04 1.671908e-04
    ## 42 1.658108e-04 1.699016e-04 1.598246e-04 1.656867e-04 1.624733e-04
    ## 43 1.691089e-04 1.734017e-04 1.612727e-04 1.684675e-04 1.648785e-04
    ## 44 1.556127e-04 1.605824e-04 1.424562e-04 1.532661e-04 1.485705e-04
    ## 45 1.596026e-04 1.642504e-04 1.491646e-04 1.581011e-04 1.539618e-04
    ## 46 1.618936e-04 1.690006e-04 1.134731e-04 1.482863e-04 1.377025e-04
    ## 47 1.854326e-04 1.922312e-04 1.366784e-04 1.712055e-04 1.608251e-04
    ## 48 1.804953e-04 1.872522e-04 1.318069e-04 1.662268e-04 1.558864e-04
    ## 49 6.165251e-05 5.569815e-05 1.000397e-04 6.939011e-05 7.865265e-05
    ## 50 8.276421e-05 8.764181e-05 4.856491e-05 6.927946e-05 6.218084e-05
    ## 51 7.000447e-05 7.491537e-05 3.829853e-05 5.655223e-05 4.969220e-05
    ## 52 5.701233e-05 5.386685e-05 8.197951e-05 5.849774e-05 6.507160e-05
    ## 53 6.043606e-05 5.767736e-05 8.268112e-05 6.096291e-05 6.693583e-05
    ## 54 5.333510e-05 5.408366e-05 5.880666e-05 4.719352e-05 4.908123e-05
    ## 55 1.893226e-05 2.615066e-05 3.244110e-05 1.126371e-05 1.168276e-05
    ## 56 6.716227e-05 6.693430e-05 7.368498e-05 6.233380e-05 6.467380e-05
    ##              26           27           28           29           30
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17                                                                 
    ## 18                                                                 
    ## 19                                                                 
    ## 20                                                                 
    ## 21                                                                 
    ## 22                                                                 
    ## 23                                                                 
    ## 24                                                                 
    ## 25                                                                 
    ## 26                                                                 
    ## 27 1.054733e-05                                                    
    ## 28 5.482471e-05 4.452986e-05                                       
    ## 29 8.598275e-05 7.553484e-05 3.151074e-05                          
    ## 30 6.314619e-05 5.263760e-05 1.087130e-05 2.313738e-05             
    ## 31 4.781174e-05 3.727701e-05 1.054880e-05 3.856310e-05 1.545524e-05
    ## 32 6.284344e-05 5.282385e-05 9.702391e-06 2.600380e-05 1.333687e-05
    ## 33 1.145827e-04 1.044652e-04 7.025446e-05 4.886901e-05 5.941584e-05
    ## 34 1.056721e-04 9.581876e-05 6.524382e-05 4.900871e-05 5.472362e-05
    ## 35 1.857959e-04 1.777255e-04 1.565832e-04 1.409896e-04 1.465886e-04
    ## 36 1.596064e-04 1.526179e-04 1.389310e-04 1.300286e-04 1.301085e-04
    ## 37 1.685867e-04 1.609540e-04 1.431131e-04 1.307261e-04 1.336206e-04
    ## 38 1.710992e-04 1.628988e-04 1.414507e-04 1.262657e-04 1.314815e-04
    ## 39 1.581675e-04 1.497913e-04 1.277371e-04 1.127870e-04 1.177613e-04
    ## 40 1.553576e-04 1.471198e-04 1.261023e-04 1.121099e-04 1.162718e-04
    ## 41 1.628257e-04 1.548424e-04 1.351053e-04 1.216290e-04 1.253937e-04
    ## 42 1.581496e-04 1.500987e-04 1.301479e-04 1.167837e-04 1.204373e-04
    ## 43 1.603916e-04 1.526228e-04 1.343216e-04 1.220459e-04 1.248178e-04
    ## 44 1.436523e-04 1.367815e-04 1.246500e-04 1.178952e-04 1.162276e-04
    ## 45 1.492589e-04 1.419125e-04 1.267149e-04 1.174012e-04 1.177505e-04
    ## 46 1.325409e-04 1.376321e-04 1.719215e-04 1.943993e-04 1.739760e-04
    ## 47 1.563553e-04 1.630259e-04 2.010010e-04 2.255972e-04 2.042231e-04
    ## 48 1.514998e-04 1.583410e-04 1.967926e-04 2.217754e-04 2.002268e-04
    ## 49 8.471825e-05 8.585285e-05 9.372305e-05 1.156160e-04 1.040198e-04
    ## 50 6.159186e-05 7.204300e-05 1.163789e-04 1.475691e-04 1.246660e-04
    ## 51 4.946387e-05 5.998716e-05 1.040924e-04 1.354002e-04 1.126098e-04
    ## 52 7.088495e-05 7.543279e-05 9.716287e-05 1.242761e-04 1.080265e-04
    ## 53 7.261882e-05 7.773523e-05 1.011154e-04 1.286026e-04 1.119862e-04
    ## 54 5.366132e-05 6.151640e-05 9.543348e-05 1.258677e-04 1.059485e-04
    ## 55 1.309496e-05 7.050823e-06 4.287671e-05 7.433327e-05 5.194926e-05
    ## 56 6.928254e-05 7.705594e-05 1.094407e-04 1.392888e-04 1.201085e-04
    ##              31           32           33           34           35
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17                                                                 
    ## 18                                                                 
    ## 19                                                                 
    ## 20                                                                 
    ## 21                                                                 
    ## 22                                                                 
    ## 23                                                                 
    ## 24                                                                 
    ## 25                                                                 
    ## 26                                                                 
    ## 27                                                                 
    ## 28                                                                 
    ## 29                                                                 
    ## 30                                                                 
    ## 31                                                                 
    ## 32 2.021296e-05                                                    
    ## 33 7.072521e-05 7.002567e-05                                       
    ## 34 6.407416e-05 6.661846e-05 1.202272e-05                          
    ## 35 1.532920e-04 1.590427e-04 9.254050e-05 9.271449e-05             
    ## 36 1.336291e-04 1.433049e-04 8.596072e-05 8.150220e-05 3.563845e-05
    ## 37 1.388570e-04 1.465112e-04 8.395053e-05 8.172367e-05 2.011302e-05
    ## 38 1.381524e-04 1.439785e-04 7.809467e-05 7.777119e-05 1.514517e-05
    ## 39 1.245101e-04 1.302681e-04 6.493299e-05 6.414436e-05 2.884647e-05
    ## 40 1.225677e-04 1.289228e-04 6.480534e-05 6.323159e-05 3.080896e-05
    ## 41 1.312906e-04 1.381351e-04 7.439371e-05 7.270272e-05 2.303317e-05
    ## 42 1.263523e-04 1.331846e-04 6.974101e-05 6.782974e-05 2.764647e-05
    ## 43 1.301199e-04 1.377085e-04 7.557633e-05 7.303721e-05 2.629335e-05
    ## 44 1.188724e-04 1.295193e-04 7.663338e-05 7.047202e-05 4.891232e-05
    ## 45 1.216700e-04 1.309072e-04 7.357113e-05 6.888879e-05 3.971744e-05
    ## 46 1.614511e-04 1.816238e-04 1.913432e-04 1.793346e-04 2.042852e-04
    ## 47 1.907559e-04 2.106463e-04 2.259536e-04 2.139310e-04 2.401633e-04
    ## 48 1.866019e-04 2.064154e-04 2.231705e-04 2.111524e-04 2.395611e-04
    ## 49 9.739917e-05 9.247363e-05 1.624794e-04 1.587266e-04 2.502448e-04
    ## 50 1.092606e-04 1.241126e-04 1.726675e-04 1.625682e-04 2.319604e-04
    ## 51 9.726376e-05 1.116626e-04 1.619675e-04 1.521665e-04 2.243826e-04
    ## 52 9.771789e-05 9.883818e-05 1.673967e-04 1.616771e-04 2.500639e-04
    ## 53 1.013359e-04 1.030532e-04 1.712913e-04 1.653731e-04 2.532614e-04
    ## 54 9.274239e-05 9.990903e-05 1.633261e-04 1.557899e-04 2.391916e-04
    ## 55 3.705608e-05 5.037359e-05 1.062922e-04 9.822029e-05 1.820662e-04
    ## 56 1.073800e-04 1.132859e-04 1.780869e-04 1.708173e-04 2.547663e-04
    ##              36           37           38           39           40
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17                                                                 
    ## 18                                                                 
    ## 19                                                                 
    ## 20                                                                 
    ## 21                                                                 
    ## 22                                                                 
    ## 23                                                                 
    ## 24                                                                 
    ## 25                                                                 
    ## 26                                                                 
    ## 27                                                                 
    ## 28                                                                 
    ## 29                                                                 
    ## 30                                                                 
    ## 31                                                                 
    ## 32                                                                 
    ## 33                                                                 
    ## 34                                                                 
    ## 35                                                                 
    ## 36                                                                 
    ## 37 1.612165e-05                                                    
    ## 38 2.895446e-05 1.370511e-05                                       
    ## 39 2.999735e-05 2.036643e-05 1.372027e-05                          
    ## 40 2.684346e-05 1.930171e-05 1.578935e-05 4.250077e-06             
    ## 41 2.129860e-05 9.960698e-06 9.696689e-06 1.043131e-05 9.589668e-06
    ## 42 2.242452e-05 1.423792e-05 1.336222e-05 7.574934e-06 5.143147e-06
    ## 43 1.620413e-05 8.803431e-06 1.458546e-05 1.413387e-05 1.175683e-05
    ## 44 1.605107e-05 2.884932e-05 3.852217e-05 3.348111e-05 2.931508e-05
    ## 45 1.262998e-05 2.009710e-05 2.844174e-05 2.354098e-05 1.948400e-05
    ## 46 1.686487e-04 1.845002e-04 1.954945e-04 1.891085e-04 1.848706e-04
    ## 47 2.045449e-04 2.204592e-04 2.315877e-04 2.252996e-04 2.210606e-04
    ## 48 2.039236e-04 2.197618e-04 2.306828e-04 2.240910e-04 2.198613e-04
    ## 49 2.308219e-04 2.362362e-04 2.351035e-04 2.214010e-04 2.196631e-04
    ## 50 2.010438e-04 2.129529e-04 2.182895e-04 2.065842e-04 2.031586e-04
    ## 51 1.943878e-04 2.057091e-04 2.104268e-04 1.983940e-04 1.950958e-04
    ## 52 2.271523e-04 2.344231e-04 2.349915e-04 2.214906e-04 2.192551e-04
    ## 53 2.299059e-04 2.374235e-04 2.382135e-04 2.247538e-04 2.224581e-04
    ## 54 2.132526e-04 2.221883e-04 2.243967e-04 2.113075e-04 2.086225e-04
    ## 55 1.577588e-04 1.656533e-04 1.671470e-04 1.539025e-04 1.513576e-04
    ## 56 2.288683e-04 2.378018e-04 2.399539e-04 2.268367e-04 2.241756e-04
    ##              41           42           43           44           45
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17                                                                 
    ## 18                                                                 
    ## 19                                                                 
    ## 20                                                                 
    ## 21                                                                 
    ## 22                                                                 
    ## 23                                                                 
    ## 24                                                                 
    ## 25                                                                 
    ## 26                                                                 
    ## 27                                                                 
    ## 28                                                                 
    ## 29                                                                 
    ## 30                                                                 
    ## 31                                                                 
    ## 32                                                                 
    ## 33                                                                 
    ## 34                                                                 
    ## 35                                                                 
    ## 36                                                                 
    ## 37                                                                 
    ## 38                                                                 
    ## 39                                                                 
    ## 40                                                                 
    ## 41                                                                 
    ## 42 4.957952e-06                                                    
    ## 43 5.351593e-06 6.717515e-06                                       
    ## 44 2.899368e-05 2.716907e-05 2.394229e-05                          
    ## 45 1.882819e-05 1.695988e-05 1.394434e-05 1.024375e-05             
    ## 46 1.859904e-04 1.838054e-04 1.809101e-04 1.569986e-04 1.671820e-04
    ## 47 2.221080e-04 2.199670e-04 2.170068e-04 1.931143e-04 2.033112e-04
    ## 48 2.211544e-04 2.188981e-04 2.160973e-04 1.921709e-04 2.023351e-04
    ## 49 2.285350e-04 2.235849e-04 2.274861e-04 2.157788e-04 2.189911e-04
    ## 50 2.091348e-04 2.049356e-04 2.056973e-04 1.851906e-04 1.928692e-04
    ## 51 2.014367e-04 1.971114e-04 1.982101e-04 1.784168e-04 1.856923e-04
    ## 52 2.275772e-04 2.227045e-04 2.258646e-04 2.115069e-04 2.159953e-04
    ## 53 2.306988e-04 2.258423e-04 2.289010e-04 2.141964e-04 2.188583e-04
    ## 54 2.162721e-04 2.115530e-04 2.139436e-04 1.972837e-04 2.029128e-04
    ## 55 1.592815e-04 1.544895e-04 1.572320e-04 1.420157e-04 1.468382e-04
    ## 56 2.318581e-04 2.271315e-04 2.295488e-04 2.128958e-04 2.185340e-04
    ##              46           47           48           49           50
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17                                                                 
    ## 18                                                                 
    ## 19                                                                 
    ## 20                                                                 
    ## 21                                                                 
    ## 22                                                                 
    ## 23                                                                 
    ## 24                                                                 
    ## 25                                                                 
    ## 26                                                                 
    ## 27                                                                 
    ## 28                                                                 
    ## 29                                                                 
    ## 30                                                                 
    ## 31                                                                 
    ## 32                                                                 
    ## 33                                                                 
    ## 34                                                                 
    ## 35                                                                 
    ## 36                                                                 
    ## 37                                                                 
    ## 38                                                                 
    ## 39                                                                 
    ## 40                                                                 
    ## 41                                                                 
    ## 42                                                                 
    ## 43                                                                 
    ## 44                                                                 
    ## 45                                                                 
    ## 46                                                                 
    ## 47 3.619665e-05                                                    
    ## 48 3.528001e-05 5.808305e-06                                       
    ## 49 2.104417e-04 2.276017e-04 2.221091e-04                          
    ## 50 1.092464e-04 1.198159e-04 1.141370e-04 1.096104e-04             
    ## 51 1.150312e-04 1.286136e-04 1.230434e-04 9.930329e-05 1.275986e-05
    ## 52 1.871065e-04 2.017943e-04 1.961918e-04 2.921575e-05 8.245215e-05
    ## 53 1.860276e-04 1.999081e-04 1.942718e-04 3.328443e-05 8.026575e-05
    ## 54 1.580704e-04 1.720115e-04 1.663996e-04 5.732146e-05 5.268679e-05
    ## 55 1.439634e-04 1.687369e-04 1.639654e-04 7.902166e-05 7.373897e-05
    ## 56 1.679219e-04 1.792460e-04 1.735320e-04 5.699968e-05 5.958053e-05
    ##              51           52           53           54           55
    ## 2                                                                  
    ## 3                                                                  
    ## 4                                                                  
    ## 5                                                                  
    ## 6                                                                  
    ## 7                                                                  
    ## 8                                                                  
    ## 9                                                                  
    ## 10                                                                 
    ## 11                                                                 
    ## 12                                                                 
    ## 13                                                                 
    ## 14                                                                 
    ## 15                                                                 
    ## 16                                                                 
    ## 17                                                                 
    ## 18                                                                 
    ## 19                                                                 
    ## 20                                                                 
    ## 21                                                                 
    ## 22                                                                 
    ## 23                                                                 
    ## 24                                                                 
    ## 25                                                                 
    ## 26                                                                 
    ## 27                                                                 
    ## 28                                                                 
    ## 29                                                                 
    ## 30                                                                 
    ## 31                                                                 
    ## 32                                                                 
    ## 33                                                                 
    ## 34                                                                 
    ## 35                                                                 
    ## 36                                                                 
    ## 37                                                                 
    ## 38                                                                 
    ## 39                                                                 
    ## 40                                                                 
    ## 41                                                                 
    ## 42                                                                 
    ## 43                                                                 
    ## 44                                                                 
    ## 45                                                                 
    ## 46                                                                 
    ## 47                                                                 
    ## 48                                                                 
    ## 49                                                                 
    ## 50                                                                 
    ## 51                                                                 
    ## 52 7.323629e-05                                                    
    ## 53 7.158781e-05 4.929567e-06                                       
    ## 54 4.355021e-05 2.980950e-05 2.805217e-05                          
    ## 55 6.131543e-05 6.964518e-05 7.219188e-05 5.782461e-05             
    ## 56 5.300567e-05 2.781226e-05 2.385183e-05 1.562159e-05 7.316100e-05

    # Convert to a full matrix if needed
    dist_matrix <- as.matrix(euclidean_distances)
    print(dist_matrix)

    ##               1            2            3            4            5
    ## 1  0.000000e+00 2.228977e-05 1.923686e-05 6.300033e-06 1.922449e-05
    ## 2  2.228977e-05 0.000000e+00 4.634166e-06 1.979748e-05 5.807367e-06
    ## 3  1.923686e-05 4.634166e-06 0.000000e+00 1.587047e-05 1.351288e-06
    ## 4  6.300033e-06 1.979748e-05 1.587047e-05 0.000000e+00 1.545900e-05
    ## 5  1.922449e-05 5.807367e-06 1.351288e-06 1.545900e-05 0.000000e+00
    ## 6  2.510122e-05 4.476402e-06 5.882784e-06 2.165374e-05 6.210366e-06
    ## 7  1.947767e-05 7.354475e-06 3.084638e-06 1.522315e-05 1.735537e-06
    ## 8  9.439574e-06 1.911386e-05 1.481852e-05 3.140692e-06 1.416699e-05
    ## 9  2.991668e-05 7.861282e-06 1.215316e-05 2.765858e-05 1.304162e-05
    ## 10 2.277153e-05 1.748207e-05 1.340338e-05 1.670609e-05 1.205229e-05
    ## 11 2.362780e-05 2.710495e-05 2.259257e-05 1.755577e-05 2.132125e-05
    ## 12 2.527583e-05 3.762699e-05 3.300297e-05 2.132719e-05 3.194465e-05
    ## 13 4.379720e-05 4.005196e-05 3.662231e-05 3.754365e-05 3.529671e-05
    ## 14 5.133393e-05 5.835646e-05 5.393444e-05 4.652001e-05 5.263378e-05
    ## 15 4.495374e-05 4.931212e-05 4.503152e-05 3.953313e-05 4.370153e-05
    ## 16 4.993929e-05 5.504723e-05 5.074683e-05 4.472129e-05 4.942026e-05
    ## 17 5.890550e-05 5.992038e-05 5.606701e-05 5.315132e-05 5.471612e-05
    ## 18 7.975462e-05 7.719937e-05 7.395145e-05 7.373090e-05 7.263601e-05
    ## 19 5.000466e-05 5.516987e-05 5.086564e-05 4.479810e-05 4.953968e-05
    ## 20 6.358186e-05 6.880090e-05 6.456305e-05 5.852670e-05 6.322866e-05
    ## 21 5.779112e-05 6.937987e-05 6.477550e-05 5.425869e-05 6.357374e-05
    ## 22 5.272705e-05 6.609410e-05 6.146043e-05 4.974139e-05 6.032889e-05
    ## 23 9.205171e-05 9.268110e-05 8.904030e-05 8.641622e-05 8.769526e-05
    ## 24 6.345088e-05 7.069919e-05 6.630927e-05 5.881536e-05 6.500099e-05
    ## 25 7.194949e-05 7.691207e-05 7.272897e-05 6.690757e-05 7.138914e-05
    ## 26 7.801185e-05 8.262707e-05 7.849629e-05 7.294185e-05 7.715241e-05
    ## 27 8.035659e-05 8.756099e-05 8.321945e-05 7.582839e-05 8.190145e-05
    ## 28 9.360660e-05 1.092390e-04 1.046115e-04 9.173190e-05 1.035327e-04
    ## 29 1.178172e-04 1.358390e-04 1.312923e-04 1.170912e-04 1.303404e-04
    ## 30 1.042354e-04 1.200887e-04 1.154627e-04 1.024812e-04 1.143886e-04
    ## 31 9.608249e-05 1.100113e-04 1.053817e-04 9.355751e-05 1.042265e-04
    ## 32 9.353299e-05 1.106091e-04 1.060222e-04 9.229639e-05 1.050223e-04
    ## 33 1.633966e-04 1.794929e-04 1.748648e-04 1.618349e-04 1.737823e-04
    ## 34 1.588487e-04 1.739156e-04 1.692816e-04 1.568704e-04 1.681511e-04
    ## 35 2.493735e-04 2.624470e-04 2.578582e-04 2.467617e-04 2.566399e-04
    ## 36 2.285705e-04 2.394862e-04 2.349961e-04 2.252857e-04 2.337204e-04
    ## 37 2.347088e-04 2.468028e-04 2.422523e-04 2.317773e-04 2.410058e-04
    ## 38 2.342328e-04 2.473729e-04 2.427804e-04 2.316340e-04 2.415655e-04
    ## 39 2.205819e-04 2.338679e-04 2.292695e-04 2.180226e-04 2.280607e-04
    ## 40 2.186471e-04 2.316381e-04 2.270492e-04 2.159890e-04 2.258310e-04
    ## 41 2.273238e-04 2.399628e-04 2.353880e-04 2.245581e-04 2.341580e-04
    ## 42 2.223938e-04 2.350900e-04 2.305124e-04 2.196430e-04 2.292846e-04
    ## 43 2.260161e-04 2.382472e-04 2.336893e-04 2.231212e-04 2.324475e-04
    ## 44 2.132664e-04 2.238215e-04 2.193496e-04 2.098682e-04 2.180671e-04
    ## 45 2.169744e-04 2.283546e-04 2.238365e-04 2.138157e-04 2.225728e-04
    ## 46 2.009674e-04 1.949436e-04 1.925562e-04 1.948147e-04 1.913418e-04
    ## 47 2.172996e-04 2.079031e-04 2.062146e-04 2.110107e-04 2.051170e-04
    ## 48 2.117678e-04 2.022267e-04 2.005655e-04 2.054760e-04 1.994732e-04
    ## 49 1.134716e-05 3.208976e-05 2.977788e-05 1.734842e-05 2.998493e-05
    ## 50 9.879039e-05 8.808968e-05 8.648808e-05 9.250326e-05 8.542075e-05
    ## 51 8.878308e-05 8.011840e-05 7.801800e-05 8.248478e-05 7.686154e-05
    ## 52 1.787987e-05 1.238559e-05 7.843749e-06 1.247902e-05 6.582884e-06
    ## 53 2.206941e-05 1.005754e-05 6.571021e-06 1.711843e-05 5.316563e-06
    ## 54 4.627957e-05 3.692360e-05 3.452570e-05 4.003693e-05 3.334339e-05
    ## 55 7.379809e-05 8.188837e-05 7.747103e-05 6.943608e-05 7.617010e-05
    ## 56 4.568302e-05 2.935771e-05 2.888854e-05 4.015700e-05 2.813694e-05
    ##               6            7            8            9           10
    ## 1  2.510122e-05 1.947767e-05 9.439574e-06 2.991668e-05 2.277153e-05
    ## 2  4.476402e-06 7.354475e-06 1.911386e-05 7.861282e-06 1.748207e-05
    ## 3  5.882784e-06 3.084638e-06 1.481852e-05 1.215316e-05 1.340338e-05
    ## 4  2.165374e-05 1.522315e-05 3.140692e-06 2.765858e-05 1.670609e-05
    ## 5  6.210366e-06 1.735537e-06 1.416699e-05 1.304162e-05 1.205229e-05
    ## 6  0.000000e+00 6.889024e-06 2.035136e-05 7.731828e-06 1.542617e-05
    ## 7  6.889024e-06 0.000000e+00 1.360905e-05 1.419886e-05 1.032568e-05
    ## 8  2.035136e-05 1.360905e-05 0.000000e+00 2.691281e-05 1.372341e-05
    ## 9  7.731828e-06 1.419886e-05 2.691281e-05 0.000000e+00 2.313776e-05
    ## 10 1.542617e-05 1.032568e-05 1.372341e-05 2.313776e-05 0.000000e+00
    ## 11 2.579653e-05 1.975868e-05 1.464167e-05 3.352128e-05 1.072606e-05
    ## 12 3.742227e-05 3.070988e-05 1.991348e-05 4.490864e-05 2.367035e-05
    ## 13 3.675963e-05 3.357498e-05 3.444919e-05 4.363029e-05 2.366969e-05
    ## 14 5.651632e-05 5.100937e-05 4.433511e-05 6.419552e-05 4.110375e-05
    ## 15 4.715918e-05 4.202502e-05 3.699721e-05 5.476207e-05 3.186771e-05
    ## 16 5.291297e-05 4.774986e-05 4.229703e-05 6.051442e-05 3.761380e-05
    ## 17 5.701753e-05 5.298294e-05 5.038032e-05 6.417512e-05 4.267786e-05
    ## 18 7.363234e-05 7.092237e-05 7.078069e-05 7.996501e-05 6.103481e-05
    ## 19 5.304370e-05 4.787043e-05 4.238065e-05 6.064776e-05 3.773975e-05
    ## 20 6.649227e-05 6.154252e-05 5.616966e-05 7.400860e-05 5.132601e-05
    ## 21 6.829219e-05 6.210599e-05 5.278329e-05 7.601548e-05 5.305580e-05
    ## 22 6.543471e-05 5.896774e-05 4.858792e-05 7.308961e-05 5.057397e-05
    ## 23 8.946450e-05 8.596012e-05 8.369341e-05 9.624695e-05 7.573512e-05
    ## 24 6.873588e-05 6.336242e-05 5.670186e-05 7.637708e-05 5.336259e-05
    ## 25 7.448013e-05 6.969174e-05 6.454528e-05 8.192701e-05 5.943020e-05
    ## 26 8.009102e-05 7.544606e-05 7.055597e-05 8.747282e-05 6.515501e-05
    ## 27 8.544449e-05 8.024514e-05 7.374617e-05 9.302602e-05 7.014805e-05
    ## 28 1.088173e-04 1.022414e-04 9.104761e-05 1.164199e-04 9.407938e-05
    ## 29 1.360201e-04 1.292264e-04 1.169585e-04 1.433794e-04 1.220203e-04
    ## 30 1.196855e-04 1.131028e-04 1.018431e-04 1.272841e-04 1.049488e-04
    ## 31 1.091626e-04 1.028237e-04 9.253182e-05 1.168634e-04 9.403058e-05
    ## 32 1.105773e-04 1.038435e-04 9.193867e-05 1.180344e-04 9.632575e-05
    ## 33 1.790266e-04 1.724812e-04 1.612449e-04 1.866498e-04 1.641672e-04
    ## 34 1.731832e-04 1.667809e-04 1.560703e-04 1.808689e-04 1.580956e-04
    ## 35 2.611333e-04 2.551390e-04 2.456060e-04 2.688625e-04 2.457389e-04
    ## 36 2.377175e-04 2.321317e-04 2.237874e-04 2.453820e-04 2.223084e-04
    ## 37 2.452798e-04 2.394626e-04 2.304603e-04 2.529904e-04 2.298551e-04
    ## 38 2.460860e-04 2.400700e-04 2.304884e-04 2.538164e-04 2.306991e-04
    ## 39 2.326265e-04 2.265744e-04 2.169007e-04 2.403581e-04 2.172535e-04
    ## 40 2.303302e-04 2.243306e-04 2.148169e-04 2.380599e-04 2.149387e-04
    ## 41 2.385671e-04 2.326397e-04 2.233283e-04 2.462914e-04 2.231568e-04
    ## 42 2.337116e-04 2.277697e-04 2.184222e-04 2.414373e-04 2.183045e-04
    ## 43 2.367617e-04 2.309112e-04 2.218251e-04 2.444772e-04 2.213398e-04
    ## 44 2.219931e-04 2.164677e-04 2.083146e-04 2.296429e-04 2.065963e-04
    ## 45 2.266935e-04 2.210030e-04 2.123855e-04 2.343820e-04 2.112703e-04
    ## 46 1.908009e-04 1.897309e-04 1.917758e-04 1.956148e-04 1.806988e-04
    ## 47 2.035331e-04 2.036382e-04 2.078839e-04 2.072450e-04 1.956021e-04
    ## 48 1.978509e-04 1.980007e-04 2.023468e-04 2.015257e-04 1.900139e-04
    ## 49 3.548724e-05 3.046411e-05 2.042784e-05 3.917446e-05 3.403200e-05
    ## 50 8.371725e-05 8.398336e-05 8.936434e-05 8.758912e-05 7.645180e-05
    ## 51 7.589603e-05 7.532015e-05 7.934948e-05 8.058404e-05 6.700269e-05
    ## 52 1.174522e-05 5.105805e-06 1.010345e-05 1.925786e-05 6.445302e-06
    ## 53 7.872767e-06 3.733758e-06 1.495067e-05 1.560433e-05 7.600988e-06
    ## 54 3.291025e-05 3.178122e-05 3.691764e-05 3.860454e-05 2.376715e-05
    ## 55 7.997838e-05 7.454347e-05 6.745706e-05 8.762856e-05 6.459089e-05
    ## 56 2.488176e-05 2.712776e-05 3.743797e-05 2.801371e-05 2.411001e-05
    ##              11           12           13           14           15
    ## 1  2.362780e-05 2.527583e-05 4.379720e-05 5.133393e-05 4.495374e-05
    ## 2  2.710495e-05 3.762699e-05 4.005196e-05 5.835646e-05 4.931212e-05
    ## 3  2.259257e-05 3.300297e-05 3.662231e-05 5.393444e-05 4.503152e-05
    ## 4  1.755577e-05 2.132719e-05 3.754365e-05 4.652001e-05 3.953313e-05
    ## 5  2.132125e-05 3.194465e-05 3.529671e-05 5.263378e-05 4.370153e-05
    ## 6  2.579653e-05 3.742227e-05 3.675963e-05 5.651632e-05 4.715918e-05
    ## 7  1.975868e-05 3.070988e-05 3.357498e-05 5.100937e-05 4.202502e-05
    ## 8  1.464167e-05 1.991348e-05 3.444919e-05 4.433511e-05 3.699721e-05
    ## 9  3.352128e-05 4.490864e-05 4.363029e-05 6.419552e-05 5.476207e-05
    ## 10 1.072606e-05 2.367035e-05 2.366969e-05 4.110375e-05 3.186771e-05
    ## 11 0.000000e+00 1.361466e-05 2.053007e-05 3.142206e-05 2.301331e-05
    ## 12 1.361466e-05 0.000000e+00 2.889268e-05 2.648400e-05 2.201844e-05
    ## 13 2.053007e-05 2.889268e-05 0.000000e+00 2.674694e-05 1.703275e-05
    ## 14 3.142206e-05 2.648400e-05 2.674694e-05 0.000000e+00 1.008288e-05
    ## 15 2.301331e-05 2.201844e-05 1.703275e-05 1.008288e-05 0.000000e+00
    ## 16 2.858688e-05 2.601524e-05 2.157815e-05 5.350903e-06 5.753830e-06
    ## 17 3.581168e-05 3.688128e-05 2.087311e-05 1.631442e-05 1.488493e-05
    ## 18 5.617537e-05 5.894198e-05 3.737897e-05 3.710839e-05 3.704343e-05
    ## 19 2.869125e-05 2.604427e-05 2.174577e-05 5.185260e-06 5.885778e-06
    ## 20 4.252094e-05 3.899791e-05 3.263002e-05 1.259314e-05 1.955683e-05
    ## 21 4.249568e-05 3.293790e-05 4.235744e-05 1.597785e-05 2.533388e-05
    ## 22 3.985247e-05 2.871909e-05 4.321991e-05 1.883869e-05 2.647245e-05
    ## 23 6.916020e-05 6.878786e-05 5.270494e-05 4.325313e-05 4.720743e-05
    ## 24 4.384042e-05 3.831696e-05 3.671865e-05 1.243707e-05 2.167175e-05
    ## 25 5.082676e-05 4.726984e-05 3.969261e-05 2.079056e-05 2.782207e-05
    ## 26 5.674415e-05 5.334221e-05 4.475237e-05 2.686050e-05 3.373139e-05
    ## 27 6.083073e-05 5.512842e-05 5.178961e-05 2.948004e-05 3.828953e-05
    ## 28 8.335429e-05 7.161460e-05 8.447799e-05 5.776269e-05 6.752485e-05
    ## 29 1.113838e-04 9.864554e-05 1.148911e-04 8.841107e-05 9.785838e-05
    ## 30 9.422401e-05 8.246201e-05 9.501751e-05 6.827261e-05 7.812457e-05
    ## 31 8.339863e-05 7.270141e-05 8.194077e-05 5.533012e-05 6.537687e-05
    ## 32 8.564492e-05 7.315504e-05 8.893604e-05 6.256051e-05 7.190507e-05
    ## 33 1.534564e-04 1.418691e-04 1.525941e-04 1.260482e-04 1.361012e-04
    ## 34 1.474467e-04 1.364121e-04 1.452195e-04 1.189735e-04 1.290564e-04
    ## 35 2.354174e-04 2.257003e-04 2.294777e-04 2.049636e-04 2.148109e-04
    ## 36 2.124035e-04 2.039694e-04 2.041716e-04 1.812047e-04 1.906751e-04
    ## 37 2.197046e-04 2.105560e-04 2.127337e-04 1.888524e-04 1.985541e-04
    ## 38 2.203567e-04 2.105851e-04 2.146144e-04 1.899726e-04 1.998432e-04
    ## 39 2.068774e-04 1.970040e-04 2.014342e-04 1.766116e-04 1.865157e-04
    ## 40 2.046103e-04 1.949087e-04 1.988283e-04 1.741921e-04 1.840578e-04
    ## 41 2.128971e-04 2.034149e-04 2.066175e-04 1.822854e-04 1.920881e-04
    ## 42 2.080307e-04 1.985090e-04 2.018621e-04 1.774572e-04 1.872745e-04
    ## 43 2.111556e-04 2.019166e-04 2.044120e-04 1.803721e-04 1.901081e-04
    ## 44 1.967603e-04 1.885422e-04 1.882607e-04 1.654982e-04 1.749057e-04
    ## 45 2.012521e-04 1.925236e-04 1.936221e-04 1.701878e-04 1.797729e-04
    ## 46 1.773441e-04 1.803375e-04 1.574542e-04 1.562681e-04 1.583251e-04
    ## 47 1.940676e-04 1.992421e-04 1.735822e-04 1.773879e-04 1.776458e-04
    ## 48 1.885756e-04 1.938841e-04 1.680773e-04 1.722303e-04 1.723428e-04
    ## 49 3.353419e-05 3.142978e-05 5.404169e-05 5.786420e-05 5.289157e-05
    ## 50 7.648007e-05 8.436355e-05 5.617455e-05 6.918273e-05 6.522365e-05
    ## 51 6.585257e-05 7.278239e-05 4.533351e-05 5.652658e-05 5.296801e-05
    ## 52 1.474882e-05 2.568458e-05 3.010361e-05 4.610970e-05 3.730958e-05
    ## 53 1.794438e-05 2.995982e-05 3.013979e-05 4.870173e-05 3.944846e-05
    ## 54 2.551836e-05 3.686228e-05 1.097786e-05 3.767936e-05 2.800946e-05
    ## 55 5.494365e-05 4.852798e-05 4.760792e-05 2.353660e-05 3.293498e-05
    ## 56 3.149395e-05 4.490828e-05 2.561799e-05 5.205335e-05 4.204360e-05
    ##              16           17           18           19           20
    ## 1  4.993929e-05 5.890550e-05 7.975462e-05 5.000466e-05 6.358186e-05
    ## 2  5.504723e-05 5.992038e-05 7.719937e-05 5.516987e-05 6.880090e-05
    ## 3  5.074683e-05 5.606701e-05 7.395145e-05 5.086564e-05 6.456305e-05
    ## 4  4.472129e-05 5.315132e-05 7.373090e-05 4.479810e-05 5.852670e-05
    ## 5  4.942026e-05 5.471612e-05 7.263601e-05 4.953968e-05 6.322866e-05
    ## 6  5.291297e-05 5.701753e-05 7.363234e-05 5.304370e-05 6.649227e-05
    ## 7  4.774986e-05 5.298294e-05 7.092237e-05 4.787043e-05 6.154252e-05
    ## 8  4.229703e-05 5.038032e-05 7.078069e-05 4.238065e-05 5.616966e-05
    ## 9  6.051442e-05 6.417512e-05 7.996501e-05 6.064776e-05 7.400860e-05
    ## 10 3.761380e-05 4.267786e-05 6.103481e-05 3.773975e-05 5.132601e-05
    ## 11 2.858688e-05 3.581168e-05 5.617537e-05 2.869125e-05 4.252094e-05
    ## 12 2.601524e-05 3.688128e-05 5.894198e-05 2.604427e-05 3.899791e-05
    ## 13 2.157815e-05 2.087311e-05 3.737897e-05 2.174577e-05 3.263002e-05
    ## 14 5.350903e-06 1.631442e-05 3.710839e-05 5.185260e-06 1.259314e-05
    ## 15 5.753830e-06 1.488493e-05 3.704343e-05 5.885778e-06 1.955683e-05
    ## 16 0.000000e+00 1.252744e-05 3.440922e-05 1.685925e-07 1.393972e-05
    ## 17 1.252744e-05 0.000000e+00 2.220514e-05 1.260327e-05 1.378528e-05
    ## 18 3.440922e-05 2.220514e-05 0.000000e+00 3.445594e-05 2.808373e-05
    ## 19 1.685925e-07 1.260327e-05 3.445594e-05 0.000000e+00 1.384101e-05
    ## 20 1.393972e-05 1.378528e-05 2.808373e-05 1.384101e-05 0.000000e+00
    ## 21 2.132602e-05 3.073452e-05 4.827086e-05 2.116124e-05 2.018727e-05
    ## 22 2.381926e-05 3.490672e-05 5.403574e-05 2.367269e-05 2.617737e-05
    ## 23 4.277400e-05 3.335740e-05 1.856063e-05 4.274799e-05 3.102677e-05
    ## 24 1.596461e-05 1.957368e-05 3.410193e-05 1.581958e-05 6.375162e-06
    ## 25 2.228094e-05 1.931541e-05 2.613684e-05 2.218840e-05 8.380884e-06
    ## 26 2.826431e-05 2.397888e-05 2.563523e-05 2.817765e-05 1.443058e-05
    ## 27 3.253762e-05 3.186170e-05 3.616773e-05 3.240985e-05 1.923669e-05
    ## 28 6.303828e-05 6.975659e-05 8.030268e-05 6.286974e-05 5.626692e-05
    ## 29 9.374852e-05 1.011195e-04 1.115753e-04 9.358105e-05 8.772618e-05
    ## 30 7.350402e-05 7.966078e-05 8.877929e-05 7.333546e-05 6.600726e-05
    ## 31 6.037200e-05 6.538714e-05 7.344369e-05 6.020511e-05 5.161856e-05
    ## 32 6.790837e-05 7.580717e-05 8.800902e-05 6.774172e-05 6.272325e-05
    ## 33 1.310524e-04 1.351717e-04 1.390402e-04 1.308862e-04 1.214338e-04
    ## 34 1.238163e-04 1.270781e-04 1.295899e-04 1.236528e-04 1.134583e-04
    ## 35 2.091801e-04 2.094245e-04 2.053919e-04 2.090299e-04 1.968490e-04
    ## 36 1.849378e-04 1.835675e-04 1.769654e-04 1.847991e-04 1.718413e-04
    ## 37 1.928642e-04 1.924070e-04 1.873041e-04 1.927187e-04 1.801703e-04
    ## 38 1.942265e-04 1.946487e-04 1.910675e-04 1.940754e-04 1.819848e-04
    ## 39 1.809214e-04 1.815964e-04 1.786033e-04 1.807690e-04 1.688145e-04
    ## 40 1.784391e-04 1.788802e-04 1.755064e-04 1.782880e-04 1.661990e-04
    ## 41 1.864364e-04 1.864958e-04 1.823618e-04 1.862876e-04 1.739962e-04
    ## 42 1.816299e-04 1.817832e-04 1.778674e-04 1.814806e-04 1.692361e-04
    ## 43 1.844299e-04 1.841588e-04 1.794701e-04 1.842833e-04 1.718214e-04
    ## 44 1.691620e-04 1.676221e-04 1.609149e-04 1.690249e-04 1.559903e-04
    ## 45 1.740544e-04 1.731574e-04 1.675020e-04 1.739123e-04 1.611488e-04
    ## 46 1.548984e-04 1.434572e-04 1.215565e-04 1.549082e-04 1.444040e-04
    ## 47 1.751158e-04 1.629137e-04 1.407376e-04 1.751560e-04 1.664778e-04
    ## 48 1.698859e-04 1.576428e-04 1.354885e-04 1.699285e-04 1.614245e-04
    ## 49 5.734194e-05 6.743346e-05 8.892661e-05 5.738156e-05 7.041433e-05
    ## 50 6.500299e-05 5.286957e-05 3.610832e-05 6.511002e-05 6.296919e-05
    ## 51 5.245854e-05 4.021478e-05 2.383208e-05 5.256031e-05 5.023003e-05
    ## 52 4.300242e-05 4.880301e-05 6.745447e-05 4.311802e-05 5.686161e-05
    ## 53 4.519813e-05 4.986454e-05 6.742547e-05 4.532541e-05 5.887138e-05
    ## 54 3.245619e-05 2.978482e-05 4.165946e-05 3.262455e-05 4.264121e-05
    ## 55 2.721984e-05 2.889138e-05 3.763713e-05 2.707655e-05 1.524817e-05
    ## 56 4.704496e-05 4.539480e-05 5.591954e-05 4.721037e-05 5.806252e-05
    ##              21           22           23           24           25
    ## 1  5.779112e-05 5.272705e-05 9.205171e-05 6.345088e-05 7.194949e-05
    ## 2  6.937987e-05 6.609410e-05 9.268110e-05 7.069919e-05 7.691207e-05
    ## 3  6.477550e-05 6.146043e-05 8.904030e-05 6.630927e-05 7.272897e-05
    ## 4  5.425869e-05 4.974139e-05 8.641622e-05 5.881536e-05 6.690757e-05
    ## 5  6.357374e-05 6.032889e-05 8.769526e-05 6.500099e-05 7.138914e-05
    ## 6  6.829219e-05 6.543471e-05 8.946450e-05 6.873588e-05 7.448013e-05
    ## 7  6.210599e-05 5.896774e-05 8.596012e-05 6.336242e-05 6.969174e-05
    ## 8  5.278329e-05 4.858792e-05 8.369341e-05 5.670186e-05 6.454528e-05
    ## 9  7.601548e-05 7.308961e-05 9.624695e-05 7.637708e-05 8.192701e-05
    ## 10 5.305580e-05 5.057397e-05 7.573512e-05 5.336259e-05 5.943020e-05
    ## 11 4.249568e-05 3.985247e-05 6.916020e-05 4.384042e-05 5.082676e-05
    ## 12 3.293790e-05 2.871909e-05 6.878786e-05 3.831696e-05 4.726984e-05
    ## 13 4.235744e-05 4.321991e-05 5.270494e-05 3.671865e-05 3.969261e-05
    ## 14 1.597785e-05 1.883869e-05 4.325313e-05 1.243707e-05 2.079056e-05
    ## 15 2.533388e-05 2.647245e-05 4.720743e-05 2.167175e-05 2.782207e-05
    ## 16 2.132602e-05 2.381926e-05 4.277400e-05 1.596461e-05 2.228094e-05
    ## 17 3.073452e-05 3.490672e-05 3.335740e-05 1.957368e-05 1.931541e-05
    ## 18 4.827086e-05 5.403574e-05 1.856063e-05 3.410193e-05 2.613684e-05
    ## 19 2.116124e-05 2.367269e-05 4.274799e-05 1.581958e-05 2.218840e-05
    ## 20 2.018727e-05 2.617737e-05 3.102677e-05 6.375162e-06 8.380884e-06
    ## 21 0.000000e+00 7.247707e-06 4.883601e-05 1.445031e-05 2.460931e-05
    ## 22 7.247707e-06 0.000000e+00 5.578037e-05 2.102631e-05 3.147213e-05
    ## 23 4.883601e-05 5.578037e-05 0.000000e+00 3.487134e-05 2.431858e-05
    ## 24 1.445031e-05 2.102631e-05 3.487134e-05 0.000000e+00 1.058510e-05
    ## 25 2.460931e-05 3.147213e-05 2.431858e-05 1.058510e-05 0.000000e+00
    ## 26 2.936840e-05 3.646100e-05 1.969953e-05 1.612595e-05 6.074456e-06
    ## 27 2.595867e-05 3.315690e-05 2.771667e-05 1.705405e-05 1.267660e-05
    ## 28 4.227844e-05 4.351506e-05 7.132425e-05 5.019024e-05 5.455946e-05
    ## 29 7.255087e-05 7.240161e-05 1.014230e-04 8.158878e-05 8.601126e-05
    ## 30 5.297386e-05 5.438098e-05 7.828636e-05 6.010480e-05 6.355754e-05
    ## 31 4.103952e-05 4.398837e-05 6.287924e-05 4.601088e-05 4.852586e-05
    ## 32 4.663827e-05 4.641975e-05 8.008030e-05 5.644744e-05 6.197053e-05
    ## 33 1.115499e-04 1.136172e-04 1.241028e-04 1.162600e-04 1.169237e-04
    ## 34 1.051085e-04 1.078225e-04 1.140630e-04 1.085968e-04 1.084237e-04
    ## 35 1.931072e-04 1.970366e-04 1.872817e-04 1.932480e-04 1.901134e-04
    ## 36 1.710316e-04 1.758493e-04 1.585039e-04 1.690179e-04 1.644791e-04
    ## 37 1.777225e-04 1.820709e-04 1.690257e-04 1.768996e-04 1.731494e-04
    ## 38 1.780194e-04 1.819139e-04 1.730733e-04 1.783085e-04 1.753334e-04
    ## 39 1.644968e-04 1.683175e-04 1.607723e-04 1.650294e-04 1.622848e-04
    ## 40 1.622996e-04 1.662563e-04 1.575900e-04 1.625198e-04 1.595648e-04
    ## 41 1.707004e-04 1.748170e-04 1.642657e-04 1.704878e-04 1.671908e-04
    ## 42 1.658108e-04 1.699016e-04 1.598246e-04 1.656867e-04 1.624733e-04
    ## 43 1.691089e-04 1.734017e-04 1.612727e-04 1.684675e-04 1.648785e-04
    ## 44 1.556127e-04 1.605824e-04 1.424562e-04 1.532661e-04 1.485705e-04
    ## 45 1.596026e-04 1.642504e-04 1.491646e-04 1.581011e-04 1.539618e-04
    ## 46 1.618936e-04 1.690006e-04 1.134731e-04 1.482863e-04 1.377025e-04
    ## 47 1.854326e-04 1.922312e-04 1.366784e-04 1.712055e-04 1.608251e-04
    ## 48 1.804953e-04 1.872522e-04 1.318069e-04 1.662268e-04 1.558864e-04
    ## 49 6.165251e-05 5.569815e-05 1.000397e-04 6.939011e-05 7.865265e-05
    ## 50 8.276421e-05 8.764181e-05 4.856491e-05 6.927946e-05 6.218084e-05
    ## 51 7.000447e-05 7.491537e-05 3.829853e-05 5.655223e-05 4.969220e-05
    ## 52 5.701233e-05 5.386685e-05 8.197951e-05 5.849774e-05 6.507160e-05
    ## 53 6.043606e-05 5.767736e-05 8.268112e-05 6.096291e-05 6.693583e-05
    ## 54 5.333510e-05 5.408366e-05 5.880666e-05 4.719352e-05 4.908123e-05
    ## 55 1.893226e-05 2.615066e-05 3.244110e-05 1.126371e-05 1.168276e-05
    ## 56 6.716227e-05 6.693430e-05 7.368498e-05 6.233380e-05 6.467380e-05
    ##              26           27           28           29           30
    ## 1  7.801185e-05 8.035659e-05 9.360660e-05 1.178172e-04 1.042354e-04
    ## 2  8.262707e-05 8.756099e-05 1.092390e-04 1.358390e-04 1.200887e-04
    ## 3  7.849629e-05 8.321945e-05 1.046115e-04 1.312923e-04 1.154627e-04
    ## 4  7.294185e-05 7.582839e-05 9.173190e-05 1.170912e-04 1.024812e-04
    ## 5  7.715241e-05 8.190145e-05 1.035327e-04 1.303404e-04 1.143886e-04
    ## 6  8.009102e-05 8.544449e-05 1.088173e-04 1.360201e-04 1.196855e-04
    ## 7  7.544606e-05 8.024514e-05 1.022414e-04 1.292264e-04 1.131028e-04
    ## 8  7.055597e-05 7.374617e-05 9.104761e-05 1.169585e-04 1.018431e-04
    ## 9  8.747282e-05 9.302602e-05 1.164199e-04 1.433794e-04 1.272841e-04
    ## 10 6.515501e-05 7.014805e-05 9.407938e-05 1.220203e-04 1.049488e-04
    ## 11 5.674415e-05 6.083073e-05 8.335429e-05 1.113838e-04 9.422401e-05
    ## 12 5.334221e-05 5.512842e-05 7.161460e-05 9.864554e-05 8.246201e-05
    ## 13 4.475237e-05 5.178961e-05 8.447799e-05 1.148911e-04 9.501751e-05
    ## 14 2.686050e-05 2.948004e-05 5.776269e-05 8.841107e-05 6.827261e-05
    ## 15 3.373139e-05 3.828953e-05 6.752485e-05 9.785838e-05 7.812457e-05
    ## 16 2.826431e-05 3.253762e-05 6.303828e-05 9.374852e-05 7.350402e-05
    ## 17 2.397888e-05 3.186170e-05 6.975659e-05 1.011195e-04 7.966078e-05
    ## 18 2.563523e-05 3.616773e-05 8.030268e-05 1.115753e-04 8.877929e-05
    ## 19 2.817765e-05 3.240985e-05 6.286974e-05 9.358105e-05 7.333546e-05
    ## 20 1.443058e-05 1.923669e-05 5.626692e-05 8.772618e-05 6.600726e-05
    ## 21 2.936840e-05 2.595867e-05 4.227844e-05 7.255087e-05 5.297386e-05
    ## 22 3.646100e-05 3.315690e-05 4.351506e-05 7.240161e-05 5.438098e-05
    ## 23 1.969953e-05 2.771667e-05 7.132425e-05 1.014230e-04 7.828636e-05
    ## 24 1.612595e-05 1.705405e-05 5.019024e-05 8.158878e-05 6.010480e-05
    ## 25 6.074456e-06 1.267660e-05 5.455946e-05 8.601126e-05 6.355754e-05
    ## 26 0.000000e+00 1.054733e-05 5.482471e-05 8.598275e-05 6.314619e-05
    ## 27 1.054733e-05 0.000000e+00 4.452986e-05 7.553484e-05 5.263760e-05
    ## 28 5.482471e-05 4.452986e-05 0.000000e+00 3.151074e-05 1.087130e-05
    ## 29 8.598275e-05 7.553484e-05 3.151074e-05 0.000000e+00 2.313738e-05
    ## 30 6.314619e-05 5.263760e-05 1.087130e-05 2.313738e-05 0.000000e+00
    ## 31 4.781174e-05 3.727701e-05 1.054880e-05 3.856310e-05 1.545524e-05
    ## 32 6.284344e-05 5.282385e-05 9.702391e-06 2.600380e-05 1.333687e-05
    ## 33 1.145827e-04 1.044652e-04 7.025446e-05 4.886901e-05 5.941584e-05
    ## 34 1.056721e-04 9.581876e-05 6.524382e-05 4.900871e-05 5.472362e-05
    ## 35 1.857959e-04 1.777255e-04 1.565832e-04 1.409896e-04 1.465886e-04
    ## 36 1.596064e-04 1.526179e-04 1.389310e-04 1.300286e-04 1.301085e-04
    ## 37 1.685867e-04 1.609540e-04 1.431131e-04 1.307261e-04 1.336206e-04
    ## 38 1.710992e-04 1.628988e-04 1.414507e-04 1.262657e-04 1.314815e-04
    ## 39 1.581675e-04 1.497913e-04 1.277371e-04 1.127870e-04 1.177613e-04
    ## 40 1.553576e-04 1.471198e-04 1.261023e-04 1.121099e-04 1.162718e-04
    ## 41 1.628257e-04 1.548424e-04 1.351053e-04 1.216290e-04 1.253937e-04
    ## 42 1.581496e-04 1.500987e-04 1.301479e-04 1.167837e-04 1.204373e-04
    ## 43 1.603916e-04 1.526228e-04 1.343216e-04 1.220459e-04 1.248178e-04
    ## 44 1.436523e-04 1.367815e-04 1.246500e-04 1.178952e-04 1.162276e-04
    ## 45 1.492589e-04 1.419125e-04 1.267149e-04 1.174012e-04 1.177505e-04
    ## 46 1.325409e-04 1.376321e-04 1.719215e-04 1.943993e-04 1.739760e-04
    ## 47 1.563553e-04 1.630259e-04 2.010010e-04 2.255972e-04 2.042231e-04
    ## 48 1.514998e-04 1.583410e-04 1.967926e-04 2.217754e-04 2.002268e-04
    ## 49 8.471825e-05 8.585285e-05 9.372305e-05 1.156160e-04 1.040198e-04
    ## 50 6.159186e-05 7.204300e-05 1.163789e-04 1.475691e-04 1.246660e-04
    ## 51 4.946387e-05 5.998716e-05 1.040924e-04 1.354002e-04 1.126098e-04
    ## 52 7.088495e-05 7.543279e-05 9.716287e-05 1.242761e-04 1.080265e-04
    ## 53 7.261882e-05 7.773523e-05 1.011154e-04 1.286026e-04 1.119862e-04
    ## 54 5.366132e-05 6.151640e-05 9.543348e-05 1.258677e-04 1.059485e-04
    ## 55 1.309496e-05 7.050823e-06 4.287671e-05 7.433327e-05 5.194926e-05
    ## 56 6.928254e-05 7.705594e-05 1.094407e-04 1.392888e-04 1.201085e-04
    ##              31           32           33           34           35
    ## 1  9.608249e-05 9.353299e-05 1.633966e-04 1.588487e-04 2.493735e-04
    ## 2  1.100113e-04 1.106091e-04 1.794929e-04 1.739156e-04 2.624470e-04
    ## 3  1.053817e-04 1.060222e-04 1.748648e-04 1.692816e-04 2.578582e-04
    ## 4  9.355751e-05 9.229639e-05 1.618349e-04 1.568704e-04 2.467617e-04
    ## 5  1.042265e-04 1.050223e-04 1.737823e-04 1.681511e-04 2.566399e-04
    ## 6  1.091626e-04 1.105773e-04 1.790266e-04 1.731832e-04 2.611333e-04
    ## 7  1.028237e-04 1.038435e-04 1.724812e-04 1.667809e-04 2.551390e-04
    ## 8  9.253182e-05 9.193867e-05 1.612449e-04 1.560703e-04 2.456060e-04
    ## 9  1.168634e-04 1.180344e-04 1.866498e-04 1.808689e-04 2.688625e-04
    ## 10 9.403058e-05 9.632575e-05 1.641672e-04 1.580956e-04 2.457389e-04
    ## 11 8.339863e-05 8.564492e-05 1.534564e-04 1.474467e-04 2.354174e-04
    ## 12 7.270141e-05 7.315504e-05 1.418691e-04 1.364121e-04 2.257003e-04
    ## 13 8.194077e-05 8.893604e-05 1.525941e-04 1.452195e-04 2.294777e-04
    ## 14 5.533012e-05 6.256051e-05 1.260482e-04 1.189735e-04 2.049636e-04
    ## 15 6.537687e-05 7.190507e-05 1.361012e-04 1.290564e-04 2.148109e-04
    ## 16 6.037200e-05 6.790837e-05 1.310524e-04 1.238163e-04 2.091801e-04
    ## 17 6.538714e-05 7.580717e-05 1.351717e-04 1.270781e-04 2.094245e-04
    ## 18 7.344369e-05 8.800902e-05 1.390402e-04 1.295899e-04 2.053919e-04
    ## 19 6.020511e-05 6.774172e-05 1.308862e-04 1.236528e-04 2.090299e-04
    ## 20 5.161856e-05 6.272325e-05 1.214338e-04 1.134583e-04 1.968490e-04
    ## 21 4.103952e-05 4.663827e-05 1.115499e-04 1.051085e-04 1.931072e-04
    ## 22 4.398837e-05 4.641975e-05 1.136172e-04 1.078225e-04 1.970366e-04
    ## 23 6.287924e-05 8.008030e-05 1.241028e-04 1.140630e-04 1.872817e-04
    ## 24 4.601088e-05 5.644744e-05 1.162600e-04 1.085968e-04 1.932480e-04
    ## 25 4.852586e-05 6.197053e-05 1.169237e-04 1.084237e-04 1.901134e-04
    ## 26 4.781174e-05 6.284344e-05 1.145827e-04 1.056721e-04 1.857959e-04
    ## 27 3.727701e-05 5.282385e-05 1.044652e-04 9.581876e-05 1.777255e-04
    ## 28 1.054880e-05 9.702391e-06 7.025446e-05 6.524382e-05 1.565832e-04
    ## 29 3.856310e-05 2.600380e-05 4.886901e-05 4.900871e-05 1.409896e-04
    ## 30 1.545524e-05 1.333687e-05 5.941584e-05 5.472362e-05 1.465886e-04
    ## 31 0.000000e+00 2.021296e-05 7.072521e-05 6.407416e-05 1.532920e-04
    ## 32 2.021296e-05 0.000000e+00 7.002567e-05 6.661846e-05 1.590427e-04
    ## 33 7.072521e-05 7.002567e-05 0.000000e+00 1.202272e-05 9.254050e-05
    ## 34 6.407416e-05 6.661846e-05 1.202272e-05 0.000000e+00 9.271449e-05
    ## 35 1.532920e-04 1.590427e-04 9.254050e-05 9.271449e-05 0.000000e+00
    ## 36 1.336291e-04 1.433049e-04 8.596072e-05 8.150220e-05 3.563845e-05
    ## 37 1.388570e-04 1.465112e-04 8.395053e-05 8.172367e-05 2.011302e-05
    ## 38 1.381524e-04 1.439785e-04 7.809467e-05 7.777119e-05 1.514517e-05
    ## 39 1.245101e-04 1.302681e-04 6.493299e-05 6.414436e-05 2.884647e-05
    ## 40 1.225677e-04 1.289228e-04 6.480534e-05 6.323159e-05 3.080896e-05
    ## 41 1.312906e-04 1.381351e-04 7.439371e-05 7.270272e-05 2.303317e-05
    ## 42 1.263523e-04 1.331846e-04 6.974101e-05 6.782974e-05 2.764647e-05
    ## 43 1.301199e-04 1.377085e-04 7.557633e-05 7.303721e-05 2.629335e-05
    ## 44 1.188724e-04 1.295193e-04 7.663338e-05 7.047202e-05 4.891232e-05
    ## 45 1.216700e-04 1.309072e-04 7.357113e-05 6.888879e-05 3.971744e-05
    ## 46 1.614511e-04 1.816238e-04 1.913432e-04 1.793346e-04 2.042852e-04
    ## 47 1.907559e-04 2.106463e-04 2.259536e-04 2.139310e-04 2.401633e-04
    ## 48 1.866019e-04 2.064154e-04 2.231705e-04 2.111524e-04 2.395611e-04
    ## 49 9.739917e-05 9.247363e-05 1.624794e-04 1.587266e-04 2.502448e-04
    ## 50 1.092606e-04 1.241126e-04 1.726675e-04 1.625682e-04 2.319604e-04
    ## 51 9.726376e-05 1.116626e-04 1.619675e-04 1.521665e-04 2.243826e-04
    ## 52 9.771789e-05 9.883818e-05 1.673967e-04 1.616771e-04 2.500639e-04
    ## 53 1.013359e-04 1.030532e-04 1.712913e-04 1.653731e-04 2.532614e-04
    ## 54 9.274239e-05 9.990903e-05 1.633261e-04 1.557899e-04 2.391916e-04
    ## 55 3.705608e-05 5.037359e-05 1.062922e-04 9.822029e-05 1.820662e-04
    ## 56 1.073800e-04 1.132859e-04 1.780869e-04 1.708173e-04 2.547663e-04
    ##              36           37           38           39           40
    ## 1  2.285705e-04 2.347088e-04 2.342328e-04 2.205819e-04 2.186471e-04
    ## 2  2.394862e-04 2.468028e-04 2.473729e-04 2.338679e-04 2.316381e-04
    ## 3  2.349961e-04 2.422523e-04 2.427804e-04 2.292695e-04 2.270492e-04
    ## 4  2.252857e-04 2.317773e-04 2.316340e-04 2.180226e-04 2.159890e-04
    ## 5  2.337204e-04 2.410058e-04 2.415655e-04 2.280607e-04 2.258310e-04
    ## 6  2.377175e-04 2.452798e-04 2.460860e-04 2.326265e-04 2.303302e-04
    ## 7  2.321317e-04 2.394626e-04 2.400700e-04 2.265744e-04 2.243306e-04
    ## 8  2.237874e-04 2.304603e-04 2.304884e-04 2.169007e-04 2.148169e-04
    ## 9  2.453820e-04 2.529904e-04 2.538164e-04 2.403581e-04 2.380599e-04
    ## 10 2.223084e-04 2.298551e-04 2.306991e-04 2.172535e-04 2.149387e-04
    ## 11 2.124035e-04 2.197046e-04 2.203567e-04 2.068774e-04 2.046103e-04
    ## 12 2.039694e-04 2.105560e-04 2.105851e-04 1.970040e-04 1.949087e-04
    ## 13 2.041716e-04 2.127337e-04 2.146144e-04 2.014342e-04 1.988283e-04
    ## 14 1.812047e-04 1.888524e-04 1.899726e-04 1.766116e-04 1.741921e-04
    ## 15 1.906751e-04 1.985541e-04 1.998432e-04 1.865157e-04 1.840578e-04
    ## 16 1.849378e-04 1.928642e-04 1.942265e-04 1.809214e-04 1.784391e-04
    ## 17 1.835675e-04 1.924070e-04 1.946487e-04 1.815964e-04 1.788802e-04
    ## 18 1.769654e-04 1.873041e-04 1.910675e-04 1.786033e-04 1.755064e-04
    ## 19 1.847991e-04 1.927187e-04 1.940754e-04 1.807690e-04 1.782880e-04
    ## 20 1.718413e-04 1.801703e-04 1.819848e-04 1.688145e-04 1.661990e-04
    ## 21 1.710316e-04 1.777225e-04 1.780194e-04 1.644968e-04 1.622996e-04
    ## 22 1.758493e-04 1.820709e-04 1.819139e-04 1.683175e-04 1.662563e-04
    ## 23 1.585039e-04 1.690257e-04 1.730733e-04 1.607723e-04 1.575900e-04
    ## 24 1.690179e-04 1.768996e-04 1.783085e-04 1.650294e-04 1.625198e-04
    ## 25 1.644791e-04 1.731494e-04 1.753334e-04 1.622848e-04 1.595648e-04
    ## 26 1.596064e-04 1.685867e-04 1.710992e-04 1.581675e-04 1.553576e-04
    ## 27 1.526179e-04 1.609540e-04 1.628988e-04 1.497913e-04 1.471198e-04
    ## 28 1.389310e-04 1.431131e-04 1.414507e-04 1.277371e-04 1.261023e-04
    ## 29 1.300286e-04 1.307261e-04 1.262657e-04 1.127870e-04 1.121099e-04
    ## 30 1.301085e-04 1.336206e-04 1.314815e-04 1.177613e-04 1.162718e-04
    ## 31 1.336291e-04 1.388570e-04 1.381524e-04 1.245101e-04 1.225677e-04
    ## 32 1.433049e-04 1.465112e-04 1.439785e-04 1.302681e-04 1.289228e-04
    ## 33 8.596072e-05 8.395053e-05 7.809467e-05 6.493299e-05 6.480534e-05
    ## 34 8.150220e-05 8.172367e-05 7.777119e-05 6.414436e-05 6.323159e-05
    ## 35 3.563845e-05 2.011302e-05 1.514517e-05 2.884647e-05 3.080896e-05
    ## 36 0.000000e+00 1.612165e-05 2.895446e-05 2.999735e-05 2.684346e-05
    ## 37 1.612165e-05 0.000000e+00 1.370511e-05 2.036643e-05 1.930171e-05
    ## 38 2.895446e-05 1.370511e-05 0.000000e+00 1.372027e-05 1.578935e-05
    ## 39 2.999735e-05 2.036643e-05 1.372027e-05 0.000000e+00 4.250077e-06
    ## 40 2.684346e-05 1.930171e-05 1.578935e-05 4.250077e-06 0.000000e+00
    ## 41 2.129860e-05 9.960698e-06 9.696689e-06 1.043131e-05 9.589668e-06
    ## 42 2.242452e-05 1.423792e-05 1.336222e-05 7.574934e-06 5.143147e-06
    ## 43 1.620413e-05 8.803431e-06 1.458546e-05 1.413387e-05 1.175683e-05
    ## 44 1.605107e-05 2.884932e-05 3.852217e-05 3.348111e-05 2.931508e-05
    ## 45 1.262998e-05 2.009710e-05 2.844174e-05 2.354098e-05 1.948400e-05
    ## 46 1.686487e-04 1.845002e-04 1.954945e-04 1.891085e-04 1.848706e-04
    ## 47 2.045449e-04 2.204592e-04 2.315877e-04 2.252996e-04 2.210606e-04
    ## 48 2.039236e-04 2.197618e-04 2.306828e-04 2.240910e-04 2.198613e-04
    ## 49 2.308219e-04 2.362362e-04 2.351035e-04 2.214010e-04 2.196631e-04
    ## 50 2.010438e-04 2.129529e-04 2.182895e-04 2.065842e-04 2.031586e-04
    ## 51 1.943878e-04 2.057091e-04 2.104268e-04 1.983940e-04 1.950958e-04
    ## 52 2.271523e-04 2.344231e-04 2.349915e-04 2.214906e-04 2.192551e-04
    ## 53 2.299059e-04 2.374235e-04 2.382135e-04 2.247538e-04 2.224581e-04
    ## 54 2.132526e-04 2.221883e-04 2.243967e-04 2.113075e-04 2.086225e-04
    ## 55 1.577588e-04 1.656533e-04 1.671470e-04 1.539025e-04 1.513576e-04
    ## 56 2.288683e-04 2.378018e-04 2.399539e-04 2.268367e-04 2.241756e-04
    ##              41           42           43           44           45
    ## 1  2.273238e-04 2.223938e-04 2.260161e-04 2.132664e-04 2.169744e-04
    ## 2  2.399628e-04 2.350900e-04 2.382472e-04 2.238215e-04 2.283546e-04
    ## 3  2.353880e-04 2.305124e-04 2.336893e-04 2.193496e-04 2.238365e-04
    ## 4  2.245581e-04 2.196430e-04 2.231212e-04 2.098682e-04 2.138157e-04
    ## 5  2.341580e-04 2.292846e-04 2.324475e-04 2.180671e-04 2.225728e-04
    ## 6  2.385671e-04 2.337116e-04 2.367617e-04 2.219931e-04 2.266935e-04
    ## 7  2.326397e-04 2.277697e-04 2.309112e-04 2.164677e-04 2.210030e-04
    ## 8  2.233283e-04 2.184222e-04 2.218251e-04 2.083146e-04 2.123855e-04
    ## 9  2.462914e-04 2.414373e-04 2.444772e-04 2.296429e-04 2.343820e-04
    ## 10 2.231568e-04 2.183045e-04 2.213398e-04 2.065963e-04 2.112703e-04
    ## 11 2.128971e-04 2.080307e-04 2.111556e-04 1.967603e-04 2.012521e-04
    ## 12 2.034149e-04 1.985090e-04 2.019166e-04 1.885422e-04 1.925236e-04
    ## 13 2.066175e-04 2.018621e-04 2.044120e-04 1.882607e-04 1.936221e-04
    ## 14 1.822854e-04 1.774572e-04 1.803721e-04 1.654982e-04 1.701878e-04
    ## 15 1.920881e-04 1.872745e-04 1.901081e-04 1.749057e-04 1.797729e-04
    ## 16 1.864364e-04 1.816299e-04 1.844299e-04 1.691620e-04 1.740544e-04
    ## 17 1.864958e-04 1.817832e-04 1.841588e-04 1.676221e-04 1.731574e-04
    ## 18 1.823618e-04 1.778674e-04 1.794701e-04 1.609149e-04 1.675020e-04
    ## 19 1.862876e-04 1.814806e-04 1.842833e-04 1.690249e-04 1.739123e-04
    ## 20 1.739962e-04 1.692361e-04 1.718214e-04 1.559903e-04 1.611488e-04
    ## 21 1.707004e-04 1.658108e-04 1.691089e-04 1.556127e-04 1.596026e-04
    ## 22 1.748170e-04 1.699016e-04 1.734017e-04 1.605824e-04 1.642504e-04
    ## 23 1.642657e-04 1.598246e-04 1.612727e-04 1.424562e-04 1.491646e-04
    ## 24 1.704878e-04 1.656867e-04 1.684675e-04 1.532661e-04 1.581011e-04
    ## 25 1.671908e-04 1.624733e-04 1.648785e-04 1.485705e-04 1.539618e-04
    ## 26 1.628257e-04 1.581496e-04 1.603916e-04 1.436523e-04 1.492589e-04
    ## 27 1.548424e-04 1.500987e-04 1.526228e-04 1.367815e-04 1.419125e-04
    ## 28 1.351053e-04 1.301479e-04 1.343216e-04 1.246500e-04 1.267149e-04
    ## 29 1.216290e-04 1.167837e-04 1.220459e-04 1.178952e-04 1.174012e-04
    ## 30 1.253937e-04 1.204373e-04 1.248178e-04 1.162276e-04 1.177505e-04
    ## 31 1.312906e-04 1.263523e-04 1.301199e-04 1.188724e-04 1.216700e-04
    ## 32 1.381351e-04 1.331846e-04 1.377085e-04 1.295193e-04 1.309072e-04
    ## 33 7.439371e-05 6.974101e-05 7.557633e-05 7.663338e-05 7.357113e-05
    ## 34 7.270272e-05 6.782974e-05 7.303721e-05 7.047202e-05 6.888879e-05
    ## 35 2.303317e-05 2.764647e-05 2.629335e-05 4.891232e-05 3.971744e-05
    ## 36 2.129860e-05 2.242452e-05 1.620413e-05 1.605107e-05 1.262998e-05
    ## 37 9.960698e-06 1.423792e-05 8.803431e-06 2.884932e-05 2.009710e-05
    ## 38 9.696689e-06 1.336222e-05 1.458546e-05 3.852217e-05 2.844174e-05
    ## 39 1.043131e-05 7.574934e-06 1.413387e-05 3.348111e-05 2.354098e-05
    ## 40 9.589668e-06 5.143147e-06 1.175683e-05 2.931508e-05 1.948400e-05
    ## 41 0.000000e+00 4.957952e-06 5.351593e-06 2.899368e-05 1.882819e-05
    ## 42 4.957952e-06 0.000000e+00 6.717515e-06 2.716907e-05 1.695988e-05
    ## 43 5.351593e-06 6.717515e-06 0.000000e+00 2.394229e-05 1.394434e-05
    ## 44 2.899368e-05 2.716907e-05 2.394229e-05 0.000000e+00 1.024375e-05
    ## 45 1.882819e-05 1.695988e-05 1.394434e-05 1.024375e-05 0.000000e+00
    ## 46 1.859904e-04 1.838054e-04 1.809101e-04 1.569986e-04 1.671820e-04
    ## 47 2.221080e-04 2.199670e-04 2.170068e-04 1.931143e-04 2.033112e-04
    ## 48 2.211544e-04 2.188981e-04 2.160973e-04 1.921709e-04 2.023351e-04
    ## 49 2.285350e-04 2.235849e-04 2.274861e-04 2.157788e-04 2.189911e-04
    ## 50 2.091348e-04 2.049356e-04 2.056973e-04 1.851906e-04 1.928692e-04
    ## 51 2.014367e-04 1.971114e-04 1.982101e-04 1.784168e-04 1.856923e-04
    ## 52 2.275772e-04 2.227045e-04 2.258646e-04 2.115069e-04 2.159953e-04
    ## 53 2.306988e-04 2.258423e-04 2.289010e-04 2.141964e-04 2.188583e-04
    ## 54 2.162721e-04 2.115530e-04 2.139436e-04 1.972837e-04 2.029128e-04
    ## 55 1.592815e-04 1.544895e-04 1.572320e-04 1.420157e-04 1.468382e-04
    ## 56 2.318581e-04 2.271315e-04 2.295488e-04 2.128958e-04 2.185340e-04
    ##              46           47           48           49           50
    ## 1  2.009674e-04 2.172996e-04 2.117678e-04 1.134716e-05 9.879039e-05
    ## 2  1.949436e-04 2.079031e-04 2.022267e-04 3.208976e-05 8.808968e-05
    ## 3  1.925562e-04 2.062146e-04 2.005655e-04 2.977788e-05 8.648808e-05
    ## 4  1.948147e-04 2.110107e-04 2.054760e-04 1.734842e-05 9.250326e-05
    ## 5  1.913418e-04 2.051170e-04 1.994732e-04 2.998493e-05 8.542075e-05
    ## 6  1.908009e-04 2.035331e-04 1.978509e-04 3.548724e-05 8.371725e-05
    ## 7  1.897309e-04 2.036382e-04 1.980007e-04 3.046411e-05 8.398336e-05
    ## 8  1.917758e-04 2.078839e-04 2.023468e-04 2.042784e-05 8.936434e-05
    ## 9  1.956148e-04 2.072450e-04 2.015257e-04 3.917446e-05 8.758912e-05
    ## 10 1.806988e-04 1.956021e-04 1.900139e-04 3.403200e-05 7.645180e-05
    ## 11 1.773441e-04 1.940676e-04 1.885756e-04 3.353419e-05 7.648007e-05
    ## 12 1.803375e-04 1.992421e-04 1.938841e-04 3.142978e-05 8.436355e-05
    ## 13 1.574542e-04 1.735822e-04 1.680773e-04 5.404169e-05 5.617455e-05
    ## 14 1.562681e-04 1.773879e-04 1.722303e-04 5.786420e-05 6.918273e-05
    ## 15 1.583251e-04 1.776458e-04 1.723428e-04 5.289157e-05 6.522365e-05
    ## 16 1.548984e-04 1.751158e-04 1.698859e-04 5.734194e-05 6.500299e-05
    ## 17 1.434572e-04 1.629137e-04 1.576428e-04 6.743346e-05 5.286957e-05
    ## 18 1.215565e-04 1.407376e-04 1.354885e-04 8.892661e-05 3.610832e-05
    ## 19 1.549082e-04 1.751560e-04 1.699285e-04 5.738156e-05 6.511002e-05
    ## 20 1.444040e-04 1.664778e-04 1.614245e-04 7.041433e-05 6.296919e-05
    ## 21 1.618936e-04 1.854326e-04 1.804953e-04 6.165251e-05 8.276421e-05
    ## 22 1.690006e-04 1.922312e-04 1.872522e-04 5.569815e-05 8.764181e-05
    ## 23 1.134731e-04 1.366784e-04 1.318069e-04 1.000397e-04 4.856491e-05
    ## 24 1.482863e-04 1.712055e-04 1.662268e-04 6.939011e-05 6.927946e-05
    ## 25 1.377025e-04 1.608251e-04 1.558864e-04 7.865265e-05 6.218084e-05
    ## 26 1.325409e-04 1.563553e-04 1.514998e-04 8.471825e-05 6.159186e-05
    ## 27 1.376321e-04 1.630259e-04 1.583410e-04 8.585285e-05 7.204300e-05
    ## 28 1.719215e-04 2.010010e-04 1.967926e-04 9.372305e-05 1.163789e-04
    ## 29 1.943993e-04 2.255972e-04 2.217754e-04 1.156160e-04 1.475691e-04
    ## 30 1.739760e-04 2.042231e-04 2.002268e-04 1.040198e-04 1.246660e-04
    ## 31 1.614511e-04 1.907559e-04 1.866019e-04 9.739917e-05 1.092606e-04
    ## 32 1.816238e-04 2.106463e-04 2.064154e-04 9.247363e-05 1.241126e-04
    ## 33 1.913432e-04 2.259536e-04 2.231705e-04 1.624794e-04 1.726675e-04
    ## 34 1.793346e-04 2.139310e-04 2.111524e-04 1.587266e-04 1.625682e-04
    ## 35 2.042852e-04 2.401633e-04 2.395611e-04 2.502448e-04 2.319604e-04
    ## 36 1.686487e-04 2.045449e-04 2.039236e-04 2.308219e-04 2.010438e-04
    ## 37 1.845002e-04 2.204592e-04 2.197618e-04 2.362362e-04 2.129529e-04
    ## 38 1.954945e-04 2.315877e-04 2.306828e-04 2.351035e-04 2.182895e-04
    ## 39 1.891085e-04 2.252996e-04 2.240910e-04 2.214010e-04 2.065842e-04
    ## 40 1.848706e-04 2.210606e-04 2.198613e-04 2.196631e-04 2.031586e-04
    ## 41 1.859904e-04 2.221080e-04 2.211544e-04 2.285350e-04 2.091348e-04
    ## 42 1.838054e-04 2.199670e-04 2.188981e-04 2.235849e-04 2.049356e-04
    ## 43 1.809101e-04 2.170068e-04 2.160973e-04 2.274861e-04 2.056973e-04
    ## 44 1.569986e-04 1.931143e-04 1.921709e-04 2.157788e-04 1.851906e-04
    ## 45 1.671820e-04 2.033112e-04 2.023351e-04 2.189911e-04 1.928692e-04
    ## 46 0.000000e+00 3.619665e-05 3.528001e-05 2.104417e-04 1.092464e-04
    ## 47 3.619665e-05 0.000000e+00 5.808305e-06 2.276017e-04 1.198159e-04
    ## 48 3.528001e-05 5.808305e-06 0.000000e+00 2.221091e-04 1.141370e-04
    ## 49 2.104417e-04 2.276017e-04 2.221091e-04 0.000000e+00 1.096104e-04
    ## 50 1.092464e-04 1.198159e-04 1.141370e-04 1.096104e-04 0.000000e+00
    ## 51 1.150312e-04 1.286136e-04 1.230434e-04 9.930329e-05 1.275986e-05
    ## 52 1.871065e-04 2.017943e-04 1.961918e-04 2.921575e-05 8.245215e-05
    ## 53 1.860276e-04 1.999081e-04 1.942718e-04 3.328443e-05 8.026575e-05
    ## 54 1.580704e-04 1.720115e-04 1.663996e-04 5.732146e-05 5.268679e-05
    ## 55 1.439634e-04 1.687369e-04 1.639654e-04 7.902166e-05 7.373897e-05
    ## 56 1.679219e-04 1.792460e-04 1.735320e-04 5.699968e-05 5.958053e-05
    ##              51           52           53           54           55
    ## 1  8.878308e-05 1.787987e-05 2.206941e-05 4.627957e-05 7.379809e-05
    ## 2  8.011840e-05 1.238559e-05 1.005754e-05 3.692360e-05 8.188837e-05
    ## 3  7.801800e-05 7.843749e-06 6.571021e-06 3.452570e-05 7.747103e-05
    ## 4  8.248478e-05 1.247902e-05 1.711843e-05 4.003693e-05 6.943608e-05
    ## 5  7.686154e-05 6.582884e-06 5.316563e-06 3.334339e-05 7.617010e-05
    ## 6  7.589603e-05 1.174522e-05 7.872767e-06 3.291025e-05 7.997838e-05
    ## 7  7.532015e-05 5.105805e-06 3.733758e-06 3.178122e-05 7.454347e-05
    ## 8  7.934948e-05 1.010345e-05 1.495067e-05 3.691764e-05 6.745706e-05
    ## 9  8.058404e-05 1.925786e-05 1.560433e-05 3.860454e-05 8.762856e-05
    ## 10 6.700269e-05 6.445302e-06 7.600988e-06 2.376715e-05 6.459089e-05
    ## 11 6.585257e-05 1.474882e-05 1.794438e-05 2.551836e-05 5.494365e-05
    ## 12 7.278239e-05 2.568458e-05 2.995982e-05 3.686228e-05 4.852798e-05
    ## 13 4.533351e-05 3.010361e-05 3.013979e-05 1.097786e-05 4.760792e-05
    ## 14 5.652658e-05 4.610970e-05 4.870173e-05 3.767936e-05 2.353660e-05
    ## 15 5.296801e-05 3.730958e-05 3.944846e-05 2.800946e-05 3.293498e-05
    ## 16 5.245854e-05 4.300242e-05 4.519813e-05 3.245619e-05 2.721984e-05
    ## 17 4.021478e-05 4.880301e-05 4.986454e-05 2.978482e-05 2.889138e-05
    ## 18 2.383208e-05 6.745447e-05 6.742547e-05 4.165946e-05 3.763713e-05
    ## 19 5.256031e-05 4.311802e-05 4.532541e-05 3.262455e-05 2.707655e-05
    ## 20 5.023003e-05 5.686161e-05 5.887138e-05 4.264121e-05 1.524817e-05
    ## 21 7.000447e-05 5.701233e-05 6.043606e-05 5.333510e-05 1.893226e-05
    ## 22 7.491537e-05 5.386685e-05 5.767736e-05 5.408366e-05 2.615066e-05
    ## 23 3.829853e-05 8.197951e-05 8.268112e-05 5.880666e-05 3.244110e-05
    ## 24 5.655223e-05 5.849774e-05 6.096291e-05 4.719352e-05 1.126371e-05
    ## 25 4.969220e-05 6.507160e-05 6.693583e-05 4.908123e-05 1.168276e-05
    ## 26 4.946387e-05 7.088495e-05 7.261882e-05 5.366132e-05 1.309496e-05
    ## 27 5.998716e-05 7.543279e-05 7.773523e-05 6.151640e-05 7.050823e-06
    ## 28 1.040924e-04 9.716287e-05 1.011154e-04 9.543348e-05 4.287671e-05
    ## 29 1.354002e-04 1.242761e-04 1.286026e-04 1.258677e-04 7.433327e-05
    ## 30 1.126098e-04 1.080265e-04 1.119862e-04 1.059485e-04 5.194926e-05
    ## 31 9.726376e-05 9.771789e-05 1.013359e-04 9.274239e-05 3.705608e-05
    ## 32 1.116626e-04 9.883818e-05 1.030532e-04 9.990903e-05 5.037359e-05
    ## 33 1.619675e-04 1.673967e-04 1.712913e-04 1.633261e-04 1.062922e-04
    ## 34 1.521665e-04 1.616771e-04 1.653731e-04 1.557899e-04 9.822029e-05
    ## 35 2.243826e-04 2.500639e-04 2.532614e-04 2.391916e-04 1.820662e-04
    ## 36 1.943878e-04 2.271523e-04 2.299059e-04 2.132526e-04 1.577588e-04
    ## 37 2.057091e-04 2.344231e-04 2.374235e-04 2.221883e-04 1.656533e-04
    ## 38 2.104268e-04 2.349915e-04 2.382135e-04 2.243967e-04 1.671470e-04
    ## 39 1.983940e-04 2.214906e-04 2.247538e-04 2.113075e-04 1.539025e-04
    ## 40 1.950958e-04 2.192551e-04 2.224581e-04 2.086225e-04 1.513576e-04
    ## 41 2.014367e-04 2.275772e-04 2.306988e-04 2.162721e-04 1.592815e-04
    ## 42 1.971114e-04 2.227045e-04 2.258423e-04 2.115530e-04 1.544895e-04
    ## 43 1.982101e-04 2.258646e-04 2.289010e-04 2.139436e-04 1.572320e-04
    ## 44 1.784168e-04 2.115069e-04 2.141964e-04 1.972837e-04 1.420157e-04
    ## 45 1.856923e-04 2.159953e-04 2.188583e-04 2.029128e-04 1.468382e-04
    ## 46 1.150312e-04 1.871065e-04 1.860276e-04 1.580704e-04 1.439634e-04
    ## 47 1.286136e-04 2.017943e-04 1.999081e-04 1.720115e-04 1.687369e-04
    ## 48 1.230434e-04 1.961918e-04 1.942718e-04 1.663996e-04 1.639654e-04
    ## 49 9.930329e-05 2.921575e-05 3.328443e-05 5.732146e-05 7.902166e-05
    ## 50 1.275986e-05 8.245215e-05 8.026575e-05 5.268679e-05 7.373897e-05
    ## 51 0.000000e+00 7.323629e-05 7.158781e-05 4.355021e-05 6.131543e-05
    ## 52 7.323629e-05 0.000000e+00 4.929567e-06 2.980950e-05 6.964518e-05
    ## 53 7.158781e-05 4.929567e-06 0.000000e+00 2.805217e-05 7.219188e-05
    ## 54 4.355021e-05 2.980950e-05 2.805217e-05 0.000000e+00 5.782461e-05
    ## 55 6.131543e-05 6.964518e-05 7.219188e-05 5.782461e-05 0.000000e+00
    ## 56 5.300567e-05 2.781226e-05 2.385183e-05 1.562159e-05 7.316100e-05
    ##              56
    ## 1  4.568302e-05
    ## 2  2.935771e-05
    ## 3  2.888854e-05
    ## 4  4.015700e-05
    ## 5  2.813694e-05
    ## 6  2.488176e-05
    ## 7  2.712776e-05
    ## 8  3.743797e-05
    ## 9  2.801371e-05
    ## 10 2.411001e-05
    ## 11 3.149395e-05
    ## 12 4.490828e-05
    ## 13 2.561799e-05
    ## 14 5.205335e-05
    ## 15 4.204360e-05
    ## 16 4.704496e-05
    ## 17 4.539480e-05
    ## 18 5.591954e-05
    ## 19 4.721037e-05
    ## 20 5.806252e-05
    ## 21 6.716227e-05
    ## 22 6.693430e-05
    ## 23 7.368498e-05
    ## 24 6.233380e-05
    ## 25 6.467380e-05
    ## 26 6.928254e-05
    ## 27 7.705594e-05
    ## 28 1.094407e-04
    ## 29 1.392888e-04
    ## 30 1.201085e-04
    ## 31 1.073800e-04
    ## 32 1.132859e-04
    ## 33 1.780869e-04
    ## 34 1.708173e-04
    ## 35 2.547663e-04
    ## 36 2.288683e-04
    ## 37 2.378018e-04
    ## 38 2.399539e-04
    ## 39 2.268367e-04
    ## 40 2.241756e-04
    ## 41 2.318581e-04
    ## 42 2.271315e-04
    ## 43 2.295488e-04
    ## 44 2.128958e-04
    ## 45 2.185340e-04
    ## 46 1.679219e-04
    ## 47 1.792460e-04
    ## 48 1.735320e-04
    ## 49 5.699968e-05
    ## 50 5.958053e-05
    ## 51 5.300567e-05
    ## 52 2.781226e-05
    ## 53 2.385183e-05
    ## 54 1.562159e-05
    ## 55 7.316100e-05
    ## 56 0.000000e+00
