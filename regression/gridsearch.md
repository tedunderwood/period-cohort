Preliminary analysis
================

``` r
knitr::opts_chunk$set(echo = TRUE, 
                      message = FALSE,
                      warning = FALSE,
                      collapse = TRUE,
                      fig.retina = 2, # Control using dpi
                      fig.width = 6,  # generated images
                      fig.pos = "t",  # pdf mode
                      fig.align = "center", 
                      out.width = "100%",
                      dev = "svg")
```

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
## [1] 1808 1848 1868 1888 1908 1928 1948 1963
## [1] 1889 1909 1929 1949 1969 1990
## [1] 1808 1832 1856 1880 1904 1928 1952 1963
## [1] 1889 1913 1937 1961 1985 1990
```

SCALE LIWC COLUMNS WITH BIZARRELY UNEVEN SCALES

``` r
for (varnum in seq(10, 89)){
  d[ , varnum] <- scale(d[ , varnum])[, 1]
}
```

``` r
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

r_squared <- function(vals, preds) {
  1 - (sum((vals - preds)^2) / sum((vals - mean(preds))^2))
}

cross_validate <- function(modelstring, data, depvar) {
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

for (varnum in seq(10, 89)){
  depvar <- colnames(d)[varnum]
  if (depvar == 'function') next  # that word breaks my as.formula!     
  bestr2 <- -0.5
  bestmodel <- 'the unknown ideal'
  bestbywidth <- 50
  bestfpwidth <- 50
  
  for (bywidth in seq(4, 24, 4)) {
    for (fpwidth in seq(4, 24, 4)) {
      bylabel = paste('by_', as.character(bywidth), sep = '')
      fplabel = paste('fp_', as.character(fpwidth), sep = '')
      modelstring <- paste(depvar, '~', 'authorage + I(authorage^2) + ', bylabel, '+', fplabel)
      thisr2 <- cross_validate(modelstring, d, depvar)
      if (thisr2 > bestr2) {
        bestr2 <- thisr2
        bestmodel <- modelstring
        bestbywidth <- bywidth
        bestfpwidth <- fpwidth
      }
    }
  }
  model <- lm(as.formula(bestmodel), data = d)
  r2 <- summary(model)$r.squared
  at <- anova_test(model, detailed = TRUE)
  cmse <- at[3,2]
  cdf <- at[3,4]
  pmse <- at[4, 2]
  pdf <- at[4,4]
  delta <- cmse / (cmse + pmse)
  adjdelta = (cmse/cdf) / ((cmse/cdf) + (pmse / pdf))
  cat(depvar, bestbywidth, bestfpwidth, delta, adjdelta, bestr2, '\n')
  varname.col <- c(varname.col, depvar)
  cmse.col <- c(cmse.col, cmse)
  pmse.col <- c(pmse.col, pmse)
  r2.col <- c(r2.col, bestr2)
  delta.col <- c(delta.col, delta)
  adjdelta.col <- c(adjdelta.col, adjdelta)
  bywidth.col <- c(bywidth.col, bestbywidth)
  fpwidth.col <- c(fpwidth.col, bestfpwidth)
  bydf.col <- c(bydf.col, cdf)
  fpdf.col <- c(fpdf.col, pdf)
}
## Analytic 24 24 0.5464306 0.4454169 0.02061517 
## Clout 20 12 0.3606309 0.4292419 0.01114944 
## Authentic 24 16 0.2733328 0.2733328 0.01922986 
## Tone 20 8 0.6338786 0.7759189 0.1269734 
## WPS 20 20 0.6326659 0.5344965 0.001060627 
## Sixltr 20 8 0.3512055 0.5198403 0.01414274 
## Dic 20 12 0.4240925 0.4954218 0.02217423 
## was_function 24 24 0.3410883 0.2565623 0.007920739 
## pronoun 24 16 0.418431 0.418431 0.0005746616 
## ppron 24 8 0.3566389 0.5257683 0.002219413 
## i 20 24 0.7511864 0.6680735 0.03591608 
## we 16 24 0.5178684 0.3494078 0.01566879 
## you 8 24 0.9645729 0.8719057 -5.680743e-05 
## shehe 24 16 0.2020363 0.2020363 0.01051649 
## they 8 8 0.6652033 0.5984202 0.01644297 
## ipron 24 20 0.9215028 0.8867011 0.02149246 
## article 12 16 0.6597894 0.5140514 0.02292833 
## prep 20 16 0.4300149 0.4300149 0.0536602 
## auxverb 12 24 0.7334331 0.5001276 0.05775332 
## adverb 20 12 0.2985899 0.3620818 0.03864422 
## conj 16 16 0.7261049 0.6653588 0.02909045 
## negate 24 12 0.5184554 0.5894125 0.02457346 
## verb 20 8 0.3132689 0.4770826 0.05658753 
## adj 16 24 0.7082268 0.5482594 0.0193073 
## compare 24 12 0.7298391 0.7827031 0.00910989 
## interrog 20 24 0.8288822 0.7635535 0.02909532 
## number 16 16 0.71632 0.6544365 0.0161284 
## quant 16 20 0.842332 0.727611 0.03345766 
## affect 12 4 0.612138 0.7749484 0.107809 
## posemo 20 16 0.6404671 0.6404671 0.1644222 
## negemo 20 16 0.5039609 0.5039609 0.01715048 
## anx 20 24 0.3501747 0.2643002 0.009528793 
## anger 24 20 0.6019672 0.5020508 0.04248378 
## sad 20 4 0.129563 0.3731952 0.09585766 
## social 24 12 0.3247008 0.3906526 0.004548516 
## family 20 20 0.6390025 0.5412987 0.01190052 
## friend 20 12 0.4751345 0.5468962 0.1103002 
## female 16 12 0.374451 0.374451 0.0109914 
## male 20 24 0.7656442 0.6853378 0.01285786 
## cogproc 20 12 0.5249973 0.5957421 0.02873581 
## insight 24 20 0.5624759 0.4615142 0.01138098 
## cause 24 20 0.6257125 0.5270739 0.05153561 
## discrep 20 12 0.4693867 0.5411754 0.02865598 
## tentat 24 24 0.5844465 0.4839027 0.03023329 
## certain 24 20 0.9233297 0.8892405 0.03218707 
## differ 16 24 0.8647384 0.7617085 0.02634147 
## percept 12 20 0.8510142 0.6750197 0.1354998 
## see 20 16 0.376777 0.376777 0.04345063 
## hear 12 20 0.7582568 0.5328389 0.1066037 
## feel 24 8 0.2952105 0.4558495 0.1068769 
## bio 12 20 0.8243174 0.6304798 0.1405375 
## body 24 24 0.8281223 0.7625865 0.1176592 
## health 24 16 0.4897252 0.4897252 0.02702824 
## sexual 24 4 0.1046219 0.3185163 0.1073786 
## ingest 8 4 0.659329 0.7437917 0.07849735 
## drives 20 24 0.5662051 0.4652857 0.03353296 
## affiliation 20 20 0.658192 0.5621229 0.03273395 
## achieve 20 4 0.2091545 0.5140625 0.01987259 
## power 24 16 0.2197681 0.2197681 0.008135988 
## reward 16 12 0.5958377 0.5958377 0.01663074 
## risk 24 16 0.5616542 0.5616542 0.04025887 
## focuspast 20 16 0.4424392 0.4424392 0.05402289 
## focuspresent 16 20 0.7855215 0.6467974 0.01442526 
## focusfuture 12 16 0.6236537 0.4747589 0.05181248 
## relativ 16 8 0.3314474 0.4264914 0.03851815 
## motion 24 20 0.4568772 0.3593041 0.03300077 
## space 20 16 0.3859761 0.3859761 0.06596926 
## time 16 8 0.560232 0.6564624 -0.0009086702 
## work 16 12 0.5877367 0.5877367 0.01462284 
## leisure 24 8 0.5983678 0.7487235 0.07192069 
## home 24 8 0.2732402 0.4292045 0.02823431 
## money 16 20 0.6029718 0.4316104 0.0109144 
## relig 20 8 0.3142078 0.4781707 0.003437434 
## death 20 12 0.3775608 0.4471402 0.01101648 
## informal 16 12 0.5964483 0.5964483 0.01931335 
## swear 20 4 0.2854006 0.6150206 0.2174132 
## netspeak 16 24 0.6751789 0.5096378 0.01384047 
## assent 16 20 0.9191875 0.8504598 0.008789228 
## nonflu 20 8 0.3140268 0.477961 0.01270677 
## filler 20 4 0.3984603 0.7259976 0.07311914

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col)
```

    ## Mean delta is  0.54207
    ## If we adjust for df it is  0.54148
    ## The average r2 is  0.04088
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.53361

<img src="gridsearch_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
