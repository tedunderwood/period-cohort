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
##  [1] 1808 1836 1840 1844 1848 1852 1856 1860 1864 1868 1872 1876 1880 1884 1888
## [16] 1892 1896 1900 1904 1908 1912 1916 1920 1924 1928 1932 1936 1940 1944 1948
## [31] 1952 1963
##  [1] 1889 1893 1897 1901 1905 1909 1913 1917 1921 1925 1929 1933 1937 1941 1945
## [16] 1949 1953 1957 1961 1965 1969 1973 1977 1981 1985 1990
##  [1] 1808 1832 1840 1848 1856 1864 1872 1880 1888 1896 1904 1912 1920 1928 1936
## [16] 1944 1952 1963
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

partial_r2s_oos <- function(modelstring, data, depvar, cohortvar, periodvar) {
  
  # See how much difference variables make out of sample.
  # Our strategy is to make predictions out of sample:
  # 1) using actual predictive variables
  # 2) with the cohort variable permuted
  # 3) with the period variable permuted
  #
  # The differences between 1 and 2 or 3 reflect lost predictive power.
  
  authors <- unique(d$author)
  authors <- sample(authors)
  tenauthsets <- chunk2(authors, 10)
  
  baseprediction <- c()
  predictionsanscohort <- c()
  predictionsansperiod <- c()
  correctyvals <- c()
  
  for (authset in tenauthsets){
    dtest <- d[d$author %in% authset,  ]
    dtrain <- d[!d$author %in% authset,  ]
    model <- lm(as.formula(modelstring), data = dtrain)
    oos_predictions <- predict(model, newdata = dtest)
    baseprediction <- c(baseprediction, oos_predictions)
    correctyvals <- c(correctyvals, dtest[depvar][[1]])
    
    cohortdata <- data.frame(dtest)
    cohortdata[cohortvar] <- sample(cohortdata[cohortvar][[1]])
    oos_predictions <- predict(model, newdata = cohortdata)
    predictionsanscohort <- c(predictionsanscohort, oos_predictions)
    
    perioddata <- data.frame(dtest)
    perioddata[periodvar] <- sample(perioddata[periodvar][[1]])
    oos_predictions <- predict(model, newdata = perioddata)
    predictionsansperiod <- c(predictionsansperiod, oos_predictions)
  }
  #cat(correctyvals, '\n')
  # cat(baseprediction, '\n')
  baser2 <- r_squared(correctyvals, baseprediction)
  sanscohortr2 <- r_squared(correctyvals, predictionsanscohort)
  sansperiodr2 <- r_squared(correctyvals, predictionsansperiod)
  
  returnvalue <- list(periodr2 = baser2 - sansperiodr2, 
                      cohortr2  = baser2 - sanscohortr2,
                      baser2 = baser2)
  returnvalue
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
pmse_oos.col <- c()
cmse_oos.col <- c()
delta_oos.col <- c()
r2_oos.col <- c()
age.col <- c()

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
  
  variances_oos <- partial_r2s_oos(bestmodel, d, depvar, bestbywidth, bestfpwidth)
  
  cmse_oos = variances_oos$cohortr2
  pmse_oos = variances_oos$periodr2
  r2_oos = variances_oos$baser2
  
  # after recording the raw mse_oos values, which may be negative
  # create non-neg values, since delta is guaranteed non-neg
  
  if (cmse_oos < 0) cmse_oos = .00000001    # debatable choice here
  if (pmse_oos < 0) pmse_oos = .00000001
  
  model <- lm(as.formula(bestmodel), data = d)
  thisr2 <- summary(model)$r.squared
  at <- anova_test(model, detailed = TRUE)
  
  cmse = 0
  pmse = 0
  cdf = 0
  pdf = 0
  
  agemse = at[1, 2] + at[2, 2]
  
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
  
  varname.col <- c(varname.col, depvar)
  delta <- cmse / (cmse + pmse)
  pmse_oos.col <- c(pmse_oos.col, pmse_oos)
  cmse_oos.col <- c(cmse_oos.col, cmse_oos)
  cmse.col <- c(cmse.col, cmse)
  pmse.col <- c(pmse.col, pmse)

  delta_oos <- cmse_oos / (cmse_oos + pmse_oos)
  
  adjdelta = (cmse/cdf) / ((cmse/cdf) + (pmse / pdf))
  cat(depvar, bestbywidth, bestfpwidth, delta, adjdelta, thisr2, '\n')
  delta_oos.col <- c(delta_oos.col, delta_oos)
  r2.col <- c(r2.col, thisr2)
  r2_oos.col <- c(r2_oos.col, r2_oos)
  delta.col <- c(delta.col, delta)
  adjdelta.col <- c(adjdelta.col, adjdelta)
  bywidth.col <- c(bywidth.col, bestbywidth)
  fpwidth.col <- c(fpwidth.col, bestfpwidth)
  bydf.col <- c(bydf.col, cdf)
  fpdf.col <- c(fpdf.col, pdf)
  age.col <- c(age.col, agemse)
}
## Analytic by_20 fp_24 0.7360433 0.6259025 0.03349767 
## Clout by_20 fp_24 0.6505936 0.5276777 0.01841032 
## Authentic by_20 fp_24 0.2596923 0.1738772 0.02851998 
## Tone by_20 fp_4 0.4383919 0.7893356 0.1440021 
## WPS by_8 fp_20 0.8031102 0.5048882 0.01464109 
## Sixltr by_16 fp_4 0.3201651 0.5855499 0.03740375 
## Dic by_20 fp_20 0.5324906 0.4767673 0.03412572 
## was_function by_20 fp_24 0.4809347 0.357295 0.0217553 
## pronoun by_24 fp_12 0.4526105 0.6231685 0.01749054 
## ppron by_24 fp_12 0.3953023 0.5666188 0.01881757 
## i by_20 fp_4 0.4220803 0.7780565 0.05087089 
## we by_20 fp_20 0.5823811 0.5273258 0.01797148 
## you by_20 fp_20 0.811188 0.7746234 0.02231925 
## shehe by_24 fp_16 0.1395728 0.1957018 0.02124582 
## they by_16 fp_24 0.8394206 0.6621957 0.02275039 
## ipron by_20 fp_12 0.8068829 0.8698786 0.03601518 
## article by_24 fp_8 0.3317048 0.5982377 0.02969788 
## prep by_20 fp_16 0.5958953 0.6389276 0.06782346 
## auxverb by_12 fp_16 0.6534719 0.5070495 0.07368519 
## adverb by_20 fp_4 0.2404928 0.6031569 0.05890832 
## conj by_16 fp_8 0.6443655 0.7310247 0.04297362 
## negate by_24 fp_16 0.3918744 0.4915071 0.03697176 
## verb by_20 fp_20 0.5504672 0.4948539 0.07326663 
## adj by_24 fp_20 0.8198453 0.8198453 0.0356062 
## compare by_24 fp_20 0.8903745 0.8903745 0.01393722 
## interrog by_20 fp_12 0.6893438 0.7802385 0.04851721 
## number by_16 fp_20 0.821626 0.6972541 0.02629761 
## quant by_24 fp_24 0.8069742 0.7581907 0.0408665 
## affect by_12 fp_24 0.9011056 0.7130589 0.1261989 
## posemo by_12 fp_16 0.8780591 0.7970635 0.1845748 
## negemo by_24 fp_20 0.4875259 0.4875259 0.02444201 
## anx by_20 fp_8 0.1485435 0.2951291 0.02371706 
## anger by_20 fp_4 0.4387467 0.789575 0.0592845 
## sad by_12 fp_24 0.5823583 0.2755146 0.1060994 
## social by_20 fp_12 0.4734472 0.589934 0.01690492 
## family by_24 fp_20 0.5536979 0.5536979 0.01838846 
## friend by_16 fp_8 0.4003085 0.5003214 0.1224197 
## female by_24 fp_12 0.2565566 0.4083486 0.02241406 
## male by_24 fp_16 0.2624558 0.3480145 0.02241734 
## cogproc by_12 fp_20 0.9268182 0.8215977 0.04526902 
## insight by_20 fp_16 0.6518388 0.691993 0.02086316 
## cause by_24 fp_12 0.4615028 0.6315456 0.05869422 
## discrep by_24 fp_24 0.5008953 0.4294487 0.0393135 
## tentat by_20 fp_8 0.391192 0.606629 0.04798186 
## certain by_12 fp_20 0.9284756 0.8251886 0.04986996 
## differ by_20 fp_20 0.3422947 0.2939599 0.03846173 
## percept by_12 fp_12 0.8284905 0.7784252 0.1538272 
## see by_20 fp_16 0.6003536 0.6431954 0.06096207 
## hear by_20 fp_24 0.8446264 0.7653497 0.112352 
## feel by_20 fp_24 0.8115054 0.7209132 0.1164225 
## bio by_24 fp_16 0.256478 0.3409889 0.1482167 
## body by_20 fp_12 0.4973265 0.6128504 0.129908 
## health by_12 fp_16 0.5538433 0.4037355 0.03653312 
## sexual by_20 fp_8 0.1585842 0.3114538 0.1173664 
## ingest by_8 fp_20 0.8498901 0.5859979 0.09806537 
## drives by_20 fp_24 0.5346023 0.4080108 0.03886858 
## affiliation by_20 fp_4 0.460543 0.8038383 0.04148737 
## achieve by_20 fp_4 0.2062604 0.5550262 0.03904032 
## power by_24 fp_8 0.1639081 0.3703257 0.019462 
## reward by_16 fp_20 0.5812881 0.4097295 0.03146656 
## risk by_12 fp_8 0.6028157 0.6234511 0.05843841 
## focuspast by_20 fp_20 0.6981761 0.6491908 0.06891901 
## focuspresent by_20 fp_20 0.8401065 0.8078154 0.03488815 
## focusfuture by_12 fp_8 0.5907511 0.6116094 0.07495936 
## relativ by_20 fp_20 0.4382149 0.3842486 0.04709023 
## motion by_16 fp_16 0.3442721 0.28252 0.04407938 
## space by_12 fp_20 0.8068649 0.6030441 0.08124868 
## time by_16 fp_24 0.6660618 0.4279053 0.0175382 
## work by_20 fp_20 0.705565 0.65719 0.02588697 
## leisure by_20 fp_16 0.1533029 0.1784909 0.07973989 
## home by_24 fp_12 0.406328 0.5778566 0.04309691 
## money by_24 fp_16 0.31447 0.4076138 0.01839019 
## relig by_24 fp_8 0.2850164 0.544606 0.01920619 
## death by_20 fp_20 0.6864327 0.6365336 0.01962599 
## informal by_16 fp_24 0.7056241 0.4733742 0.03000628 
## swear by_24 fp_8 0.1366429 0.3219455 0.2160155 
## netspeak by_16 fp_8 0.5484058 0.6455869 0.02448763 
## assent by_16 fp_8 0.5847278 0.6786724 0.02618648 
## nonflu by_16 fp_4 0.3300683 0.5964599 0.03636504 
## filler by_8 fp_20 0.7935023 0.4899696 0.09035007

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col, pmse_oos = pmse_oos.col, cmse_oos = cmse_oos.col, delta_oos = delta_oos.col, r2_oos = r2_oos.col)
write.csv(res.df, file = 'gridsearch_delta_oos2.csv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.54598
    ## If we adjust for df it is  0.5677
    ## The average r2 is  0.05395
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.55041
    ## 
    ## Mean delta measured oos is  0.48207
    ## but as a weighted average it is: 0.52328

<img src="gridsearch_oos_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
