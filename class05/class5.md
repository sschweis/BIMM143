Class 05
================

``` r
# Class 5 R graphics

# 2A. Line plot
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)


plot( weight$Age, weight$Weight, pch=15, cex=1.5, 
      typ= "b", lwd=2, ylim=c(2,10), 
      xlab="Age (months)", ylab = "Weight (kg)", 
      main = "Baby Weight with Age" )
```

![](class5_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# 2B. Barplot
feat <- read.table("bimm143_05_rstats/feature_counts.txt", 
                     sep= "\t", header=TRUE)
par(mar= c(3.1, 11.1, 4.1, 2))
barplot(height = feat$Count, horiz = TRUE, names.arg = feat$Feature,
        las=1, xlim = c(0, 80000), 
        main = "Number of features in the mouse GRCm38 genome")
```

![](class5_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# 3A. Providing Color Vectors

MF <- read.table("bimm143_05_rstats/male_female_counts.txt", sep = "\t", header= TRUE)
par(mar= c(8.1, 3.1, 4.1, 2))
barplot(height = MF$Count, names.arg = MF$Sample, las=2, col = cm.colors(10))
```

![](class5_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
