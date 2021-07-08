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
## [1] 1889 1905 1921 1937 1953 1969 1985 1990
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
  fiveauthsets <- chunk2(authors, 5)
  rsquaredvals <- c()
  
  for (authset in fiveauthsets){
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
## Analytic by_20 fp_16 0.7577409 0.7896231 0.03289917 
## Clout by_20 fp_20 0.4587939 0.4041163 0.02147811 
## Authentic by_24 fp_12 0.05596308 0.1059944 0.02838165 
## Tone by_12 fp_16 0.9593651 0.9279428 0.14804 
## WPS by_20 fp_16 0.6008864 0.6437049 0.009899119 
## Sixltr by_20 fp_8 0.4639579 0.6750356 0.03234831 
## Dic by_24 fp_20 0.2984343 0.2984343 0.03060867 
## was_function by_20 fp_20 0.5890841 0.5342057 0.01999042 
## pronoun by_24 fp_24 0.4844705 0.4134261 0.01643846 
## ppron by_24 fp_20 0.5212397 0.5212397 0.01412352 
## i by_20 fp_4 0.4220803 0.7780565 0.05087089 
## we by_20 fp_16 0.4044321 0.4490005 0.02068344 
## you by_20 fp_8 0.7688872 0.8886977 0.0230294 
## shehe by_24 fp_12 0.1237852 0.2203005 0.02330205 
## they by_20 fp_12 0.2280711 0.3209888 0.02289 
## ipron by_20 fp_12 0.8068829 0.8698786 0.03601518 
## article by_24 fp_16 0.6982861 0.7763662 0.02446178 
## prep by_20 fp_12 0.6528121 0.7505274 0.06660313 
## auxverb by_12 fp_8 0.5857077 0.606652 0.07657784 
## adverb by_20 fp_8 0.2662042 0.4654307 0.05704687 
## conj by_20 fp_20 0.8212386 0.7861072 0.03537081 
## negate by_20 fp_8 0.3670651 0.5819152 0.04188987 
## verb by_20 fp_12 0.545111 0.6572222 0.07310222 
## adj by_24 fp_8 0.6806649 0.8647647 0.03735482 
## compare by_24 fp_24 0.8535758 0.8138532 0.01447291 
## interrog by_20 fp_12 0.6893438 0.7802385 0.04851721 
## number by_16 fp_8 0.669196 0.7521321 0.0318026 
## quant by_20 fp_12 0.4435272 0.5604884 0.04761438 
## affect by_20 fp_4 0.2368163 0.5983036 0.1257349 
## posemo by_12 fp_24 0.9598524 0.8670282 0.1827736 
## negemo by_20 fp_24 0.5605093 0.4334973 0.02505866 
## anx by_20 fp_24 0.4126036 0.2964964 0.01785945 
## anger by_24 fp_4 0.2433963 0.6587235 0.05868036 
## sad by_16 fp_8 0.1437381 0.2011507 0.1087547 
## social by_8 fp_20 0.6263081 0.2952785 0.01929796 
## family by_24 fp_24 0.6553254 0.5877929 0.01765435 
## friend by_16 fp_24 0.8272043 0.642243 0.1174003 
## female by_20 fp_4 0.2648671 0.6336234 0.02713703 
## male by_20 fp_12 0.6398328 0.7397448 0.02117151 
## cogproc by_24 fp_8 0.4020238 0.6685364 0.04127732 
## insight by_20 fp_12 0.563976 0.6742164 0.02218919 
## cause by_24 fp_12 0.4615028 0.6315456 0.05869422 
## discrep by_24 fp_24 0.5008953 0.4294487 0.0393135 
## tentat by_24 fp_24 0.3517581 0.2892556 0.04469576 
## certain by_12 fp_16 0.9409065 0.8967466 0.04953282 
## differ by_12 fp_8 0.6408843 0.6606552 0.04709173 
## percept by_12 fp_20 0.8510142 0.6750197 0.1532687 
## see by_20 fp_8 0.5783697 0.7670188 0.06181422 
## hear by_20 fp_12 0.672359 0.7665402 0.1174244 
## feel by_20 fp_24 0.8115054 0.7209132 0.1164225 
## bio by_16 fp_20 0.8473821 0.7351804 0.1493041 
## body by_20 fp_24 0.8555278 0.7803669 0.1261985 
## health by_12 fp_16 0.5505197 0.4005041 0.03668952 
## sexual by_24 fp_4 0.05476225 0.2579452 0.1179274 
## ingest by_24 fp_20 0.5516243 0.5516243 0.0864737 
## drives by_24 fp_4 0.1946263 0.5918297 0.04398648 
## affiliation by_20 fp_20 0.6756517 0.6249744 0.0372291 
## achieve by_24 fp_4 0.1444841 0.503306 0.03767525 
## power by_24 fp_24 0.2420625 0.1932409 0.01424526 
## reward by_20 fp_20 0.3871636 0.3357271 0.02579302 
## risk by_12 fp_20 0.6769035 0.4324105 0.05649822 
## focuspast by_20 fp_24 0.6446023 0.52113 0.07031646 
## focuspresent by_20 fp_4 0.633962 0.8926277 0.04170528 
## focusfuture by_12 fp_12 0.6198644 0.5425267 0.07331195 
## relativ by_24 fp_8 0.08584196 0.2197913 0.04732707 
## motion by_20 fp_12 0.5188197 0.6330487 0.04530756 
## space by_20 fp_16 0.4112298 0.4559737 0.07995179 
## time by_16 fp_12 0.5967797 0.5967797 0.01928622 
## work by_20 fp_8 0.5449647 0.7418899 0.02971034 
## leisure by_24 fp_24 0.9296918 0.9084024 0.07873603 
## home by_24 fp_20 0.3703577 0.3703577 0.04212497 
## money by_20 fp_12 0.1922471 0.2757841 0.01935557 
## relig by_20 fp_20 0.2311624 0.1938941 0.01150706 
## death by_24 fp_20 0.4470037 0.4470037 0.01592768 
## informal by_16 fp_16 0.7095242 0.6468891 0.03021254 
## swear by_16 fp_12 0.4096581 0.4096581 0.2170371 
## netspeak by_24 fp_20 0.5766119 0.5766119 0.01726123 
## assent by_16 fp_8 0.5847278 0.6786724 0.02618648 
## nonflu by_20 fp_24 0.6204515 0.49516 0.02993087 
## filler by_20 fp_12 0.4386118 0.5555709 0.08131035

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col)
write.csv(res.df, file = 'crossvalidated_delta_knit.tsv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.52887
    ## If we adjust for df it is  0.57424
    ## The average r2 is  0.05334
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.52252

<img src="gridsearch_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
