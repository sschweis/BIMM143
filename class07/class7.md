Class 7: R Functions and Packages
================
Sandy Schweis
4/23/2019

\#Functions Revisted

We will source a file from online with our functions from last day

``` r
source("http://tinyurl.com/rescale-R")
```

Try out the last day’s rescale()
    function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Try the rescale2() function that catches string inputs

\#Find missing NA values into vectors

Start with a simple example of the larger problem I am trying to solve

``` r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

Try putting these together with an AND

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take the sum() to find out how many TRUE values we have and thus how
many NAs we have in both x and y

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

Now I can this into a function…

``` r
both_na <- function(x, y) {
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(x, y)
```

    ## [1] 1

``` r
both_na(x, c(NA,3, NA, 2, NA))
```

    ## [1] 2

Try both\_na2() function on our different length input vectors

``` r
both_na2(x, y)
```

    ## [1] 1

\#\#Intersect function

``` r
df1
```

    ##     IDs exp
    ## 1 gene1   2
    ## 2 gene2   1
    ## 3 gene3   1

``` r
df2
```

    ##     IDs exp
    ## 1 gene2  -2
    ## 2 gene4  NA
    ## 3 gene3   1
    ## 4 gene5   2

Make things simple

``` r
x <- df1$IDs
y <- df2$IDs

x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
intersect(x, y)
```

    ## [1] "gene2" "gene3"

``` r
which(x %in% y)
```

    ## [1] 2 3

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
y[y %in% x]
```

    ## [1] "gene2" "gene3"

``` r
cbind(x[x %in% y],
      y[y %in% x])
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

Use RStudio shortcut CODE \> EXTRACT FUNCTION to turn our snippet into a
working function

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y],
        y[y %in% x])
}
```

``` r
gene_intersect2(df1, df2)
```

    ##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
    ## 2 gene2   1                               -2
    ## 3 gene3   1                                1

``` r
gene_intersect3(df1, df2)
```

    ##     IDs exp exp2
    ## 2 gene2   1   -2
    ## 3 gene3   1    1

``` r
merge(df1, df2, by = "IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

Find average grade dropping the worst homework score

``` r
a <- c(100, 100, 100, 100, 100, 100, 100, 90)
b <- c(100, 90, 90, 90, 90, 90, 97, 80)
```

``` r
grade <- function(a) {
  (sum(a) - min(a)) / (length(a)-1)
}
```

``` r
grade(a)
```

    ## [1] 100
