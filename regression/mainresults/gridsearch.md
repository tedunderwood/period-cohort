Finding the best models through a grid search
================

Needed packages.

``` r
library(tidyverse)
library(rstatix)
```

Bring in the data. Limiting to 1890-1989. Not requiring multiple works,
but requiring us\_national. Once this chunk is done executing, there
should be 5,573 books in d.

``` r
d <- read_tsv("../liwc/liwc_w_meta.tsv")

d <- d %>%
  filter(firstpub < 1990 & firstpub > 1889) %>% 
  filter(us_national == TRUE) %>%
  rename(author = hathi_author)
```

Now we convert the numeric columns for birthyear and firstpub into
factors. Six different factors are created for each numeric variable,
using different binwidths.

To avoid discarding data from the tails, we use a function that creates
ad-hoc breaks and requires each new factor level to have at least 20
books. So the breaks toward the end, and especially at the beginning of
birthyear, may be more widely spaced than in the middle of the timeline.

``` r
discretize <- function(numericvar, width) {
  minval = min(numericvar) - 1
  maxval = max(numericvar) + 1
  breaks = c()
  for (i in seq(minval, maxval, width)) {
    breaks <- c(breaks, i)
  }
  breaks <- c(breaks, maxval)
  
  neededbreaks = c(breaks[1])
  previous = breaks[1]
  for (idx in seq(2, length(breaks) - 1)) {
    prevct = sum(numericvar >= previous & numericvar < breaks[idx])
    nextct = sum(numericvar >= breaks[idx] & numericvar < breaks[idx + 1])
    if (prevct > 20 & nextct > 20 & breaks[idx] != 1989 ){
      previous = breaks[idx]
      neededbreaks <- c(neededbreaks, breaks[idx])
    }
  }
breaks <- c(neededbreaks, maxval)

# The algorithm above ensures that breaks are less tightly packed in tails of 
# distributions. But when the list of breaks gets quite short it can be
# tricky to arrange the breaks to compromise between even chronological
# spacing and (roughly!) even numbers of volumes in factor levels.
# The following cases are therefore addressed manually:

if (minval == 1889 & width == 24) breaks = c(1889, 1914, 1938, 1961, 1990)
if (minval == 1808 & width == 24) breaks = c(1808, 1854, 1878, 1902, 1926, 1990)
if (minval == 1808 & width == 20) breaks = c(1808, 1850, 1870, 1890, 1910, 1930, 1990)
if (minval == 1889 & width == 16) breaks = c(1889, 1905, 1921, 1937, 1953, 1968, 1982, 1990)

print(breaks)
result <- cut(numericvar, breaks = as.integer(breaks), labels = as.character(breaks[1: length(breaks) -1]))
result
}

for (width in seq(4, 24, 4)){
  bylabel = paste('by_', as.character(width), sep = '')
  fplabel = paste('fp_', as.character(width), sep = '')
  d[bylabel] <- discretize(as.integer(d$birthyear), width)
  d[fplabel] <- discretize(as.integer(d$firstpub), width)
}
##  [1] 1808 1836 1840 1844 1848 1852 1856 1860 1864 1868 1872 1876 1880 1884
## [15] 1888 1892 1896 1900 1904 1908 1912 1916 1920 1924 1928 1932 1936 1940
## [29] 1944 1948 1952 1963
##  [1] 1889 1893 1897 1901 1905 1909 1913 1917 1921 1925 1929 1933 1937 1941
## [15] 1945 1949 1953 1957 1961 1965 1969 1973 1977 1981 1985 1990
##  [1] 1808 1832 1840 1848 1856 1864 1872 1880 1888 1896 1904 1912 1920 1928
## [15] 1936 1944 1952 1963
##  [1] 1889 1897 1905 1913 1921 1929 1937 1945 1953 1961 1969 1977 1985 1990
##  [1] 1808 1832 1844 1856 1868 1880 1892 1904 1916 1928 1940 1952 1963
##  [1] 1889 1901 1913 1925 1937 1949 1961 1973 1985 1990
##  [1] 1808 1840 1856 1872 1888 1904 1920 1936 1952 1963
## [1] 1889 1905 1921 1937 1953 1968 1982 1990
## [1] 1808 1850 1870 1890 1910 1930 1990
## [1] 1889 1909 1929 1949 1969 1990
## [1] 1808 1854 1878 1902 1926 1990
## [1] 1889 1914 1938 1961 1990
```

SCALE LIWC COLUMNS WITH BIZARRELY UNEVEN SCALES

All the dependent variables are converted to zscores, so they will have
comparable scales.

``` r
for (varnum in seq(10, 89)){
  d[ , varnum] <- scale(d[ , varnum])[, 1]
}
```

Now we actually do a grid search to find the best binwdiths for firstpub
and birthyear. We select the model with highest overall r2, and in doing
that use fivefold cross-validation *on unseen authors*.

After selecting the best model we find out how variance is apportioned
across variables.

``` r
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

r_squared <- function(vals, preds) {
  1 - (sum((vals - preds)^2) / sum((vals - mean(preds))^2))
}

cross_validate <- function(modelstring, data, depvar) {
  
  # Testing a model specified by modelstring.
  # We test it on out-of-sample authors.
  
  authors <- unique(d$author)
  authors <- sample(authors)
  tenauthsets <- chunk2(authors, 10)
  rsquaredvals <- c()
  
  for (authset in tenauthsets){
    dtest <- d[d$author %in% authset,  ]
    dtrain <- d[!d$author %in% authset,  ]
    model <- lm(as.formula(modelstring), data = dtrain)
    oos_predictions <- predict(model, newdata = dtest)
    r2 <- r_squared(dtest[depvar], oos_predictions)
    rsquaredvals <- c(rsquaredvals, r2)
  }
  mean(rsquaredvals)
}

periodnames = c('fp_4', 'fp_8', 'fp_12', 'fp_16', 'fp_20', 'fp_24')
cohortnames = c('by_4', 'by_8', 'by_12', 'by_16', 'by_20', 'by_24')

varname.col <- c()
cmse.col <- c()
pmse.col <- c()
r2.col <- c()
delta.col <- c()
adjdelta.col <- c()
bywidth.col <- c()
fpwidth.col <- c()
bydf.col <- c()
fpdf.col <- c()

# Iterate across dependent variables:

for (varnum in seq(10, 89)){
  depvar <- colnames(d)[varnum]
  if (depvar == 'function') next  # that word breaks my as.formula!     
  bestr2 <- -0.5
  bestmodel <- 'the unknown ideal'
  bestbywidth <- 50
  bestfpwidth <- 50
  
  # For each dependent variables, do a grid search to find the
  # granularity of factors that produce the highest r2.
  # This is assessed using the cross_validate function above,
  # which tests on out-of-sample authors.
  
  for (bywidth in seq(4, 24, 4)) {
    for (fpwidth in seq(4, 24, 4)) {
      
      bylabel = paste('by_', as.character(bywidth), sep = '')
      fplabel = paste('fp_', as.character(fpwidth), sep = '')
      modelstring <- paste(depvar, '~', 'authorage + I(authorage^2) + ', bylabel, '+', fplabel)
      thisr2 <- cross_validate(modelstring, d, depvar)
      
      if (thisr2 > bestr2) {
        bestr2 <- thisr2
        bestmodel <- modelstring
        bestbywidth <- bylabel
        bestfpwidth <- fplabel
      }
    }
  }
  
  # We're done finding the best model specification.
  # Now we train it on all data and examine sums of squares.
  
  model <- lm(as.formula(bestmodel), data = d)
  thisr2 <- summary(model)$r.squared
  at <- anova_test(model, detailed = TRUE)
  
  cmse = 0
  pmse = 0
  cdf = 0
  pdf = 0
  
  for (rownum in seq(3, dim(at)[1])) {
    variable_name = at[rownum, 1]
    if (variable_name %in% cohortnames){
      cmse = cmse + at[rownum, 2]
      cdf = cdf + at[rownum, 4]
    }
    else if (variable_name %in% periodnames) {
      pmse = pmse + at[rownum, 2]
      pdf = pdf + at[rownum, 4]
    }
    else {
      print('Model error.')
    }
  }

  delta <- cmse / (cmse + pmse)
  adjdelta = (cmse/cdf) / ((cmse/cdf) + (pmse / pdf))
  cat(depvar, bestbywidth, bestfpwidth, delta, adjdelta, thisr2, '\n')
  varname.col <- c(varname.col, depvar)
  cmse.col <- c(cmse.col, cmse)
  pmse.col <- c(pmse.col, pmse)
  r2.col <- c(r2.col, thisr2)
  delta.col <- c(delta.col, delta)
  adjdelta.col <- c(adjdelta.col, adjdelta)
  bywidth.col <- c(bywidth.col, bestbywidth)
  fpwidth.col <- c(fpwidth.col, bestfpwidth)
  bydf.col <- c(bydf.col, cdf)
  fpdf.col <- c(fpdf.col, pdf)
}
## Analytic by_20 fp_24 0.7360433 0.6259025 0.03349767 
## Clout by_20 fp_12 0.4592051 0.5760211 0.02217937 
## Authentic by_20 fp_16 0.1936713 0.2237392 0.02904634 
## Tone by_12 fp_16 0.8818388 0.8027899 0.1497175 
## WPS by_8 fp_20 0.8031102 0.5048882 0.01464109 
## Sixltr by_20 fp_20 0.5769607 0.5217776 0.02781237 
## Dic by_20 fp_24 0.6427873 0.5191549 0.0325034 
## was_function by_24 fp_16 0.3791723 0.4781145 0.02357588 
## pronoun by_24 fp_16 0.5123554 0.6118029 0.0166409 
## ppron by_24 fp_24 0.4684925 0.3979823 0.01618675 
## i by_12 fp_20 0.7729949 0.5532224 0.04846038 
## we by_24 fp_24 0.2988276 0.2422159 0.01888461 
## you by_20 fp_8 0.7688872 0.8886977 0.0230294 
## shehe by_24 fp_4 0.1371486 0.4881478 0.02853965 
## they by_16 fp_12 0.5157476 0.5157476 0.02708019 
## ipron by_20 fp_24 0.9310704 0.8901646 0.03427086 
## article by_20 fp_24 0.1979 0.1289475 0.0239536 
## prep by_20 fp_12 0.6528121 0.7505274 0.06660313 
## auxverb by_24 fp_20 0.6711497 0.6711497 0.06525669 
## adverb by_20 fp_4 0.2404928 0.6031569 0.05890832 
## conj by_20 fp_8 0.5613817 0.754404 0.03956937 
## negate by_24 fp_16 0.3918744 0.4915071 0.03697176 
## verb by_20 fp_12 0.545111 0.6572222 0.07310222 
## adj by_24 fp_20 0.8198453 0.8198453 0.0356062 
## compare by_20 fp_12 0.5786297 0.68722 0.01794668 
## interrog by_20 fp_12 0.6893438 0.7802385 0.04851721 
## number by_8 fp_16 0.7406405 0.5171111 0.03467187 
## quant by_20 fp_24 0.9517765 0.922131 0.04331421 
## affect by_20 fp_4 0.2368163 0.5983036 0.1257349 
## posemo by_12 fp_12 0.9457561 0.9269016 0.1829115 
## negemo by_20 fp_8 0.3332195 0.5453276 0.0308346 
## anx by_24 fp_20 0.3498025 0.3498025 0.01883659 
## anger by_20 fp_8 0.5778579 0.7666435 0.0563979 
## sad by_24 fp_8 0.06889894 0.1816639 0.1075913 
## social by_20 fp_20 0.52375 0.4680257 0.01565063 
## family by_24 fp_16 0.5322721 0.6305864 0.01958677 
## friend by_20 fp_8 0.2372235 0.4273936 0.1208702 
## female by_20 fp_16 0.4109245 0.4556609 0.02205012 
## male by_12 fp_20 0.8327161 0.6441448 0.02770723 
## cogproc by_24 fp_24 0.7019095 0.638469 0.03837287 
## insight by_20 fp_20 0.6532041 0.6010902 0.02091716 
## cause by_24 fp_16 0.4793865 0.5800465 0.05800575 
## discrep by_24 fp_16 0.4696002 0.570457 0.03827338 
## tentat by_24 fp_24 0.3517581 0.2892556 0.04469576 
## certain by_12 fp_12 0.8914662 0.8566024 0.05033337 
## differ by_24 fp_24 0.6346494 0.5657503 0.03680398 
## percept by_12 fp_24 0.9732085 0.9083149 0.1501584 
## see by_20 fp_16 0.6003536 0.6431954 0.06096207 
## hear by_20 fp_24 0.8446264 0.7653497 0.112352 
## feel by_20 fp_4 0.3821068 0.748005 0.1250434 
## bio by_24 fp_16 0.256478 0.3409889 0.1482167 
## body by_20 fp_16 0.6768416 0.7153713 0.127259 
## health by_12 fp_16 0.5538433 0.4037355 0.03653312 
## sexual by_20 fp_24 0.4116729 0.2956957 0.108525 
## ingest by_20 fp_4 0.2728889 0.6430442 0.09161994 
## drives by_20 fp_16 0.3874002 0.4314514 0.04170832 
## affiliation by_20 fp_16 0.5969011 0.639891 0.0382777 
## achieve by_20 fp_4 0.2062604 0.5550262 0.03904032 
## power by_20 fp_8 0.2539216 0.4495877 0.02031653 
## reward by_20 fp_8 0.2530927 0.4485042 0.03049276 
## risk by_20 fp_24 0.6142634 0.4886129 0.04258512 
## focuspast by_20 fp_24 0.6446023 0.52113 0.07031646 
## focuspresent by_20 fp_8 0.6680208 0.8284548 0.0403642 
## focusfuture by_12 fp_24 0.7068189 0.3966845 0.06968307 
## relativ by_20 fp_12 0.3841063 0.4994621 0.04847586 
## motion by_20 fp_16 0.5120543 0.5573831 0.04556267 
## space by_20 fp_16 0.4148838 0.4597149 0.07992106 
## time by_16 fp_24 0.6660618 0.4279053 0.0175382 
## work by_20 fp_16 0.6670056 0.706199 0.0259764 
## leisure by_20 fp_16 0.1533029 0.1784909 0.07973989 
## home by_24 fp_12 0.406328 0.5778566 0.04309691 
## money by_24 fp_12 0.2402978 0.387484 0.02053417 
## relig by_24 fp_8 0.2850164 0.544606 0.01920619 
## death by_24 fp_4 0.04751786 0.230373 0.02636742 
## informal by_16 fp_8 0.5466771 0.6439887 0.03669946 
## swear by_24 fp_12 0.1854572 0.3128872 0.2130136 
## netspeak by_16 fp_24 0.6008382 0.3608053 0.0232177 
## assent by_16 fp_12 0.7439668 0.7439668 0.02217122 
## nonflu by_20 fp_24 0.6204515 0.49516 0.02993087 
## filler by_20 fp_16 0.5290108 0.5740746 0.07960889

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col)
write.csv(res.df, file = 'crossvalidated_delta_knit.tsv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.52071
    ## If we adjust for df it is  0.55804
    ## The average r2 is  0.05338
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.53617

<img src="gridsearch_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
