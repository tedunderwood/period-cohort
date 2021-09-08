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
## Analytic by_20 fp_12 0.6681105 0.7630827 0.03422989 
## Clout by_20 fp_16 0.4681484 0.5136822 0.020758 
## Authentic by_24 fp_16 0.03502906 0.05163916 0.02783461 
## Tone by_20 fp_12 0.6900123 0.7807736 0.1388706 
## WPS by_16 fp_20 0.6891869 0.5257705 0.008492994 
## Sixltr by_20 fp_4 0.3690147 0.7373363 0.03755357 
## Dic by_24 fp_20 0.2984343 0.2984343 0.03060867 
## was_function by_24 fp_8 0.3663573 0.6343064 0.02543075 
## pronoun by_24 fp_12 0.4526105 0.6231685 0.01749054 
## ppron by_24 fp_24 0.4684925 0.3979823 0.01618675 
## i by_20 fp_20 0.7123425 0.6645517 0.04319582 
## we by_24 fp_24 0.2988276 0.2422159 0.01888461 
## you by_16 fp_20 0.8051269 0.6738179 0.01730634 
## shehe by_20 fp_24 0.1720505 0.1108597 0.02058208 
## they by_24 fp_12 0.1132943 0.2035298 0.02178499 
## ipron by_20 fp_24 0.9310704 0.8901646 0.03427086 
## article by_24 fp_24 0.4906389 0.4194259 0.02604751 
## prep by_20 fp_20 0.6977723 0.6487544 0.06660152 
## auxverb by_24 fp_24 0.7270885 0.6664603 0.06544314 
## adverb by_20 fp_4 0.2404928 0.6031569 0.05890832 
## conj by_16 fp_16 0.7384146 0.6791923 0.04107408 
## negate by_16 fp_20 0.6131057 0.4420709 0.04008897 
## verb by_20 fp_12 0.545111 0.6572222 0.07310222 
## adj by_24 fp_20 0.8198453 0.8198453 0.0356062 
## compare by_20 fp_20 0.964007 0.95541 0.01465408 
## interrog by_20 fp_20 0.8389646 0.806496 0.0457911 
## number by_16 fp_24 0.821529 0.6331866 0.02676859 
## quant by_24 fp_12 0.210241 0.3474367 0.04544104 
## affect by_20 fp_20 0.3910129 0.3393482 0.1168099 
## posemo by_12 fp_16 0.8780591 0.7970635 0.1845748 
## negemo by_20 fp_8 0.3332195 0.5453276 0.0308346 
## anx by_24 fp_20 0.3498025 0.3498025 0.01883659 
## anger by_24 fp_4 0.2433963 0.6587235 0.05868036 
## sad by_16 fp_20 0.2806962 0.1632615 0.1055801 
## social by_20 fp_24 0.6267778 0.5018981 0.01433184 
## family by_12 fp_24 0.8243035 0.5613143 0.02231934 
## friend by_16 fp_20 0.8786891 0.7836267 0.1171498 
## female by_20 fp_16 0.4109245 0.4556609 0.02205012 
## male by_24 fp_24 0.4164686 0.3486521 0.01882394 
## cogproc by_20 fp_12 0.580998 0.6893057 0.04158337 
## insight by_20 fp_16 0.6518388 0.691993 0.02086316 
## cause by_12 fp_8 0.6147104 0.6351019 0.0670271 
## discrep by_24 fp_24 0.5008953 0.4294487 0.0393135 
## tentat by_24 fp_24 0.3517581 0.2892556 0.04469576 
## certain by_12 fp_24 0.9642832 0.880427 0.04920978 
## differ by_24 fp_20 0.3732902 0.3732902 0.03957034 
## percept by_12 fp_20 0.8510142 0.6750197 0.1532687 
## see by_20 fp_16 0.6003536 0.6431954 0.06096207 
## hear by_20 fp_20 0.6947589 0.6455006 0.1171965 
## feel by_20 fp_20 0.6827298 0.632557 0.1173614 
## bio by_20 fp_4 0.3949512 0.7580593 0.1556478 
## body by_20 fp_8 0.4752095 0.684866 0.1298476 
## health by_24 fp_24 0.5691309 0.497656 0.02829954 
## sexual by_16 fp_8 0.1347199 0.1893269 0.1166511 
## ingest by_24 fp_20 0.5516243 0.5516243 0.0864737 
## drives by_20 fp_8 0.3789795 0.5942557 0.04187661 
## affiliation by_20 fp_16 0.5969011 0.639891 0.0382777 
## achieve by_20 fp_12 0.3105623 0.418852 0.03321257 
## power by_24 fp_8 0.1639081 0.3703257 0.019462 
## reward by_16 fp_8 0.4338714 0.5347917 0.03420547 
## risk by_24 fp_8 0.3061129 0.5696093 0.04903352 
## focuspast by_20 fp_12 0.5712613 0.6807027 0.07119626 
## focuspresent by_20 fp_20 0.8401065 0.8078154 0.03488815 
## focusfuture by_12 fp_12 0.6198644 0.5425267 0.07331195 
## relativ by_20 fp_24 0.5126842 0.3869673 0.04604293 
## motion by_20 fp_8 0.4664168 0.6771999 0.04635576 
## space by_12 fp_8 0.6377342 0.6575859 0.08327692 
## time by_16 fp_20 0.8131786 0.6851736 0.01493225 
## work by_20 fp_20 0.705565 0.65719 0.02588697 
## leisure by_12 fp_16 0.7402179 0.6084891 0.08589663 
## home by_24 fp_20 0.3703577 0.3703577 0.04212497 
## money by_12 fp_12 0.7017648 0.6311753 0.02669222 
## relig by_24 fp_20 0.3699222 0.3699222 0.01593569 
## death by_20 fp_16 0.4213131 0.4662853 0.02263308 
## informal by_16 fp_8 0.5466771 0.6439887 0.03669946 
## swear by_20 fp_16 0.4782721 0.5238208 0.2145503 
## netspeak by_16 fp_20 0.7352912 0.5813917 0.02139014 
## assent by_16 fp_16 0.7550735 0.6980807 0.02185216 
## nonflu by_20 fp_24 0.6204515 0.49516 0.02993087 
## filler by_20 fp_24 0.5778458 0.4509361 0.07932108

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col, pmse_oos = pmse_oos.col, cmse_oos = cmse_oos.col, delta_oos = delta_oos.col, r2_oos = r2_oos.col)
write.csv(res.df, file = 'gridsearch_delta_oos5.csv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.54432
    ## If we adjust for df it is  0.5578
    ## The average r2 is  0.053
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.54788
    ## 
    ## Mean delta measured oos is  0.47156
    ## but as a weighted average it is: 0.52736

<img src="gridsearch_oos_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
