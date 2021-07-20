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
## Analytic by_20 fp_16 0.764659 0.795876 0.03284901 
## Clout by_20 fp_20 0.4587939 0.4041163 0.02147811 
## Authentic by_24 fp_16 0.03502906 0.05163916 0.02783461 
## Tone by_12 fp_16 0.8818388 0.8027899 0.1497175 
## WPS by_20 fp_20 0.6397778 0.5869222 0.008909021 
## Sixltr by_20 fp_16 0.5751431 0.6189722 0.02837416 
## Dic by_24 fp_20 0.2984343 0.2984343 0.03060867 
## was_function by_24 fp_24 0.3787216 0.3137469 0.02290114 
## pronoun by_24 fp_20 0.5772421 0.5772421 0.01325376 
## ppron by_24 fp_16 0.4645593 0.5654876 0.01780649 
## i by_24 fp_24 0.507201 0.43564 0.03852332 
## we by_16 fp_24 0.4748174 0.2531949 0.02131948 
## you by_20 fp_24 0.9588364 0.9332264 0.0187032 
## shehe by_24 fp_16 0.1395728 0.1957018 0.02124582 
## they by_16 fp_16 0.6685265 0.6020098 0.02538467 
## ipron by_20 fp_24 0.9310704 0.8901646 0.03427086 
## article by_24 fp_12 0.4485489 0.6193079 0.0271834 
## prep by_20 fp_20 0.6977723 0.6487544 0.06660152 
## auxverb by_12 fp_12 0.5696806 0.4905248 0.07611998 
## adverb by_20 fp_8 0.2662042 0.4654307 0.05704687 
## conj by_20 fp_16 0.6805583 0.7188289 0.037353 
## negate by_16 fp_20 0.6131057 0.4420709 0.04008897 
## verb by_20 fp_16 0.6039494 0.6466326 0.07139238 
## adj by_24 fp_8 0.6806649 0.8647647 0.03735482 
## compare by_24 fp_12 0.5099192 0.6754258 0.0172356 
## interrog by_20 fp_12 0.6893438 0.7802385 0.04851721 
## number by_16 fp_24 0.821529 0.6331866 0.02676859 
## quant by_20 fp_12 0.4435272 0.5604884 0.04761438 
## affect by_20 fp_16 0.3717756 0.4152544 0.1182834 
## posemo by_12 fp_24 0.9598524 0.8670282 0.1827736 
## negemo by_20 fp_8 0.3332195 0.5453276 0.0308346 
## anx by_24 fp_8 0.1046065 0.2595239 0.02327088 
## anger by_20 fp_20 0.8660058 0.8379362 0.05272764 
## sad by_16 fp_4 0.1298982 0.3093314 0.1117531 
## social by_20 fp_24 0.6267778 0.5018981 0.01433184 
## family by_24 fp_12 0.4708824 0.640272 0.02001143 
## friend by_16 fp_24 0.8272043 0.642243 0.1174003 
## female by_24 fp_24 0.3846457 0.3191767 0.01939698 
## male by_20 fp_20 0.5549713 0.4994085 0.02186292 
## cogproc by_20 fp_20 0.8154887 0.7795307 0.03959351 
## insight by_20 fp_8 0.5725529 0.7627369 0.02225971 
## cause by_24 fp_4 0.3434169 0.7583502 0.06317949 
## discrep by_20 fp_12 0.2439753 0.3405142 0.0395604 
## tentat by_20 fp_24 0.5801007 0.4532275 0.0459286 
## certain by_20 fp_20 0.9489176 0.9369522 0.04025405 
## differ by_20 fp_12 0.3219841 0.4317622 0.0390881 
## percept by_12 fp_12 0.8284905 0.7784252 0.1538272 
## see by_20 fp_20 0.740616 0.6955146 0.05788423 
## hear by_20 fp_20 0.6947589 0.6455006 0.1171965 
## feel by_20 fp_8 0.4538499 0.6660429 0.1222431 
## bio by_20 fp_8 0.4843827 0.692744 0.1529773 
## body by_24 fp_12 0.09047277 0.1659331 0.1257911 
## health by_24 fp_16 0.2503989 0.333806 0.02990671 
## sexual by_20 fp_16 0.2341329 0.2683916 0.1129825 
## ingest by_8 fp_20 0.8498901 0.5859979 0.09806537 
## drives by_20 fp_20 0.58784 0.5329272 0.03707503 
## affiliation by_12 fp_20 0.7753235 0.5565119 0.04016353 
## achieve by_20 fp_4 0.2062604 0.5550262 0.03904032 
## power by_24 fp_8 0.1639081 0.3703257 0.019462 
## reward by_24 fp_20 0.2895897 0.2895897 0.02494999 
## risk by_24 fp_20 0.4061171 0.4061171 0.04721931 
## focuspast by_20 fp_8 0.5023188 0.7078043 0.07310489 
## focuspresent by_20 fp_8 0.6680208 0.8284548 0.0403642 
## focusfuture by_12 fp_4 0.5401421 0.7193164 0.07769749 
## relativ by_20 fp_20 0.4382149 0.3842486 0.04709023 
## motion by_20 fp_24 0.7162176 0.6022744 0.04318783 
## space by_20 fp_16 0.4148838 0.4597149 0.07992106 
## time by_16 fp_12 0.5967797 0.5967797 0.01928622 
## work by_20 fp_20 0.705565 0.65719 0.02588697 
## leisure by_24 fp_20 0.8288262 0.8288262 0.07923321 
## home by_24 fp_20 0.3703577 0.3703577 0.04212497 
## money by_24 fp_20 0.06993177 0.06993177 0.01862626 
## relig by_24 fp_20 0.3699222 0.3699222 0.01593569 
## death by_20 fp_16 0.4213131 0.4662853 0.02263308 
## informal by_16 fp_12 0.5964483 0.5964483 0.03424159 
## swear by_20 fp_16 0.4782721 0.5238208 0.2145503 
## netspeak by_24 fp_16 0.2976801 0.3886704 0.01873525 
## assent by_24 fp_8 0.1653922 0.3728453 0.01686944 
## nonflu by_24 fp_16 0.488455 0.5888654 0.02535901 
## filler by_20 fp_20 0.4444608 0.3902597 0.08112926

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col, pmse_oos = pmse_oos.col, cmse_oos = cmse_oos.col, delta_oos = delta_oos.col, r2_oos = r2_oos.col)
write.csv(res.df, file = 'gridsearch_delta_oos4.csv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.51756
    ## If we adjust for df it is  0.54545
    ## The average r2 is  0.0528
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.53088
    ## 
    ## Mean delta measured oos is  0.43673
    ## but as a weighted average it is: 0.51278

<img src="gridsearch_oos_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
