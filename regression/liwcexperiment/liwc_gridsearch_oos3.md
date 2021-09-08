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
## Clout by_24 fp_20 0.2792976 0.2792976 0.0195134 
## Authentic by_24 fp_16 0.03502906 0.05163916 0.02783461 
## Tone by_20 fp_16 0.6598714 0.6995262 0.1395673 
## WPS by_20 fp_16 0.5564329 0.6008527 0.01042196 
## Sixltr by_20 fp_4 0.3690147 0.7373363 0.03755357 
## Dic by_20 fp_12 0.5118069 0.6265021 0.03501975 
## was_function by_24 fp_24 0.3787216 0.3137469 0.02290114 
## pronoun by_24 fp_24 0.4844705 0.4134261 0.01643846 
## ppron by_24 fp_8 0.3986614 0.6654254 0.02005114 
## i by_20 fp_4 0.4220803 0.7780565 0.05087089 
## we by_20 fp_16 0.3970308 0.441388 0.02076364 
## you by_20 fp_20 0.811188 0.7746234 0.02231925 
## shehe by_20 fp_24 0.1720505 0.1108597 0.02058208 
## they by_16 fp_12 0.5157476 0.5157476 0.02708019 
## ipron by_20 fp_24 0.9310704 0.8901646 0.03427086 
## article by_24 fp_12 0.4485489 0.6193079 0.0271834 
## prep by_20 fp_4 0.5030292 0.8293083 0.06970726 
## auxverb by_24 fp_24 0.7270885 0.6664603 0.06544314 
## adverb by_20 fp_12 0.2971941 0.4035508 0.05458224 
## conj by_16 fp_24 0.8756338 0.7252965 0.03852652 
## negate by_20 fp_4 0.306807 0.6799466 0.04739957 
## verb by_20 fp_20 0.5504672 0.4948539 0.07326663 
## adj by_24 fp_20 0.8198453 0.8198453 0.0356062 
## compare by_24 fp_24 0.8535758 0.8138532 0.01447291 
## interrog by_20 fp_16 0.7300701 0.7644619 0.04779582 
## number by_16 fp_16 0.6826698 0.6173667 0.03129289 
## quant by_20 fp_4 0.2747265 0.6451627 0.05098715 
## affect by_16 fp_24 0.6581653 0.4192872 0.117633 
## posemo by_12 fp_20 0.9844041 0.9582506 0.1822324 
## negemo by_24 fp_16 0.5065251 0.6062478 0.02515182 
## anx by_20 fp_20 0.4064345 0.3539163 0.01963762 
## anger by_20 fp_16 0.6967797 0.7338669 0.05526479 
## sad by_20 fp_20 0.1905174 0.1584515 0.1050804 
## social by_20 fp_8 0.4510459 0.6635206 0.01762451 
## family by_24 fp_12 0.4708824 0.640272 0.02001143 
## friend by_16 fp_20 0.8786891 0.7836267 0.1171498 
## female by_24 fp_4 0.212607 0.6183327 0.02715741 
## male by_20 fp_16 0.4736728 0.5192193 0.02346211 
## cogproc by_24 fp_24 0.7019095 0.638469 0.03837287 
## insight by_20 fp_24 0.8525951 0.7763074 0.01888291 
## cause by_12 fp_12 0.6526521 0.5774376 0.06588847 
## discrep by_24 fp_24 0.5008953 0.4294487 0.0393135 
## tentat by_24 fp_8 0.2580123 0.5105701 0.04760829 
## certain by_20 fp_12 0.7727475 0.8447356 0.04151485 
## differ by_24 fp_24 0.6346494 0.5657503 0.03680398 
## percept by_12 fp_20 0.8510142 0.6750197 0.1532687 
## see by_20 fp_8 0.5783697 0.7670188 0.06181422 
## hear by_20 fp_8 0.5812795 0.7691463 0.1208138 
## feel by_20 fp_8 0.4538499 0.6660429 0.1222431 
## bio by_20 fp_20 0.7807431 0.7401712 0.1493126 
## body by_20 fp_20 0.8585448 0.8292204 0.1259515 
## health by_8 fp_16 0.7061688 0.4740283 0.03861934 
## sexual by_20 fp_16 0.2341329 0.2683916 0.1129825 
## ingest by_24 fp_16 0.6297622 0.7184249 0.08502647 
## drives by_20 fp_12 0.4179673 0.5346643 0.04093283 
## affiliation by_24 fp_16 0.220694 0.2981419 0.03555591 
## achieve by_24 fp_8 0.1871539 0.4085414 0.03333093 
## power by_24 fp_12 0.1750615 0.2979615 0.01725041 
## reward by_16 fp_16 0.6079519 0.5376857 0.02934744 
## risk by_24 fp_8 0.3061129 0.5696093 0.04903352 
## focuspast by_20 fp_20 0.6981761 0.6491908 0.06891901 
## focuspresent by_20 fp_20 0.8401065 0.8078154 0.03488815 
## focusfuture by_12 fp_24 0.7068189 0.3966845 0.06968307 
## relativ by_20 fp_24 0.5126842 0.3869673 0.04604293 
## motion by_20 fp_16 0.5120543 0.5573831 0.04556267 
## space by_20 fp_24 0.54825 0.4213526 0.07654873 
## time by_24 fp_16 0.3282387 0.4229446 0.01103225 
## work by_20 fp_4 0.4489057 0.7963316 0.03352521 
## leisure by_16 fp_12 0.81999 0.81999 0.08254871 
## home by_20 fp_20 0.1937902 0.1612832 0.03775206 
## money by_24 fp_12 0.2402978 0.387484 0.02053417 
## relig by_20 fp_8 0.1512359 0.2995436 0.01545806 
## death by_16 fp_8 0.4495768 0.5505973 0.0243968 
## informal by_16 fp_24 0.7056241 0.4733742 0.03000628 
## swear by_20 fp_8 0.2990817 0.5059482 0.2183817 
## netspeak by_24 fp_24 0.2751588 0.2216139 0.01829562 
## assent by_16 fp_16 0.7550735 0.6980807 0.02185216 
## nonflu by_20 fp_12 0.5201297 0.6342669 0.03124104 
## filler by_20 fp_20 0.4444608 0.3902597 0.08112926

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col, pmse_oos = pmse_oos.col, cmse_oos = cmse_oos.col, delta_oos = delta_oos.col, r2_oos = r2_oos.col)
write.csv(res.df, file = 'gridsearch_delta_oos3.csv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.52623
    ## If we adjust for df it is  0.57133
    ## The average r2 is  0.05295
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.53082
    ## 
    ## Mean delta measured oos is  0.47164
    ## but as a weighted average it is: 0.51435

<img src="gridsearch_oos_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
